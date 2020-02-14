#include <stdio.h>
#include <sys/mman.h>
#include <signal.h>
#include <stdbool.h>
#include <rtnet.h>
#include <math.h>
#include <rtai_sem.h>
#include <rtai_shm.h>
#include <rtai_mbx.h>

#include <semaphore.h>

#include "./network/socket_UDP.h"
#include "./utils/utils.h"
#include "./header/mv_typ.h"
#include "./header/const.h"
#include "./header/type.h"
#include "./kinematics/kinematicsR.h"
#include "./control/controlR.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h> 
#include <unistd.h>
#include <string.h>
//include da Comunicazione
#include <time.h>
//include per thread
#include <pthread.h>

#define __CONTROL_GAIN_TRANS_DISPLACEMENT__   0.2	//proporzionale traslazione
#define __CONTROL_GAIN_ROT_DISPLACEMENT__      0.2	//proporzionale rotazione
#define __TRAN_MAX__ 5.0E-4 //norma massima traslazione [m] --> 25 cm/s
#define __ROT_MAX__	 3.6E-3 //norma massima rotazione ee [rad] --> 10 °/s
#define _USE_MATH_DEFINES

// Thread info structure
typedef struct
{
    pthread_t thread_ptr;
    unsigned long thread_num;
    unsigned int thread_active;
    unsigned long shrdmem_num;
    unsigned long shrdmem_mutex_num;
} thread_param;



volatile bool stop_mod5, abnormal_behaviour;
static float af_PosErr[GBM_NAX_OPEN];    // following error [delta giri]
static float af_PosErrMax[GBM_NAX_OPEN]; // massimo following error consentito [delta giri]

// utility funtions **********************************************************************

// trasposta di matrice
int matrix_transpose(orientation_matrix matrix, orientation_matrix result) {
	memset(&(result[0][0]), 0, 9 * sizeof(double));
	int i, j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			result[i][j] = matrix[j][i];
		}
	}
	return 0;
}
//norma 2 di un vettore
double vector_norm(translation_vector v){
		
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	
}
// prodotto di matrici di rotzione
int matrix_matrix_mult(orientation_matrix matrix1, orientation_matrix matrix2, orientation_matrix result) {

	memset(&(result[0][0]), 0, 9 * sizeof(double));
	int i, j, k;
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) for (k = 0; k<3; k++) result[i][j] += matrix1[i][k] * matrix2[k][j];
	return 0;

}
// prodotto matrice di rotazione per vettore
int matrix_vector_mult(orientation_matrix matrix, translation_vector vector, translation_vector result) {

	memset(&(result[0]), 0, 3 * sizeof(double));
	int i, j;
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) result[i] += matrix[i][j] * vector[j];
	return 0;

}
//prodotto matrice omogenea per vettore traslazione
int pose_vector_mult(pose_matrix matrix, translation_vector vector, translation_vector result){
	
	memset(&(result[0]), 0, 3 * sizeof(double));
	orientation_matrix R;
	translation_vector t;
	pose_mat_getRT(matrix, R, t);
	
	matrix_vector_mult(R, vector, result);
	result[0] += t[0];
	result[1] += t[1];
	result[2] += t[2];
		
	return 0;
}


//prodotto tra due matrici omogenee
int pose_pose_mult(const pose_matrix m1, const pose_matrix m2, pose_matrix res) {
	memset(&(res[0]), 0, 4 * 3 * sizeof(double));

	double m1_full[4][4] = { {m1[0][0], m1[0][1], m1[0][2], m1[0][3]},{m1[1][0], m1[1][1], m1[1][2], m1[1][3]},{m1[2][0], m1[2][1], m1[2][2], m1[2][3]},{0,0,0,1} };
	double m2_full[4][4] = { {m2[0][0], m2[0][1], m2[0][2], m2[0][3]},{m2[1][0], m2[1][1], m2[1][2], m2[1][3]},{m2[2][0], m2[2][1], m2[2][2], m2[2][3]},{0,0,0,1} };

	int i, j, k;

	for (i = 0; i < 3; i++) for (j = 0; j < 4; j++) for (k = 0; k < 4; k++) res[i][j] += m1_full[i][k] * m2_full[k][j];

	return 0;
}
//print translation vector
void print_vector(const translation_vector v) {
	printf("\n t = (%6.2f, %6.2f, %6.2f)\n", v[0], v[1], v[2]);
}
// differenza (v1-v2) di vettori posizione
int vector_diff(const translation_vector v1, const translation_vector v2, translation_vector res) {
	memset(&(res[0]), 0, 3 * sizeof(double));
	int i;
	for (i = 0;i < 3;i++)	res[i] = v1[i] - v2[i];

	return 0;
}
//somma (v1+v2) di vettori posizione
int vector_sum(const translation_vector v1, const translation_vector v2, translation_vector res) {
	memset(&(res[0]), 0, 3 * sizeof(double));
	int i;
	for (i = 0;i < 3;i++)	res[i] = v1[i] + v2[i];

	return 0;
}
// prodotto di uno scalare per un vettore
int scalar_vect_mult(double c, const translation_vector v, translation_vector res) {
	memset(&(res[0]), 0, 3 * sizeof(double));
	int i;
	for (i = 0;i < 3;i++)	res[i] = c * v[i];

	return 0;
}
// prodotto di vettore posizione: v1' * v2
double vect_vect_mult(const translation_vector v1, const translation_vector v2) {

	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];

}
// calcolo matrice omogenea a partire dai parametri DH
void evalDHmatrix(pose_matrix A, double alpha, double a, double theta, double d, double offset) {
	memset(&(A[0][0]), 0, 12 * sizeof(double));

	double cth = cos(theta + offset);
	double sth = sin(theta + offset);
	double cal = cos(alpha);
	double sal = sin(alpha);

	A[0][0] = cth;
	A[0][1] = -sth * cal;
	A[0][2] = sth * sal;
	A[0][3] = a * cth;
	A[1][0] = sth;
	A[1][1] = cth * cal;
	A[1][2] = -cth * sal;
	A[1][3] = a * sth;
	A[2][0] = 0;
	A[2][1] = sal;
	A[2][2] = cal;
	A[2][3] = d;

}
//cinematica diretta -> ritorna tutte le posizioni dei giunti nello spazio
void fkine(joint_vector const q, translation_vector joints[]) {
	memset(&(joints[0]), 0, 3 * 6 * sizeof(double));
	joint_vector q_DH;
	pose_matrix A01, A02, A03, A04, A05, A06, A12, A23, A34, A45, A56;
	double offset_tool = 0.19; //offset lunghezza tool lungo z ee [m]

	//calcolo coord giunto con convenzione comau
	q_DH[0] = -1 * q[0];
	q_DH[1] = q[1] - 90;
	q_DH[2] = -1 * q[2] - 90;
	q_DH[3] = -1 * q[3];
	q_DH[4] = q[4];
	q_DH[5] = -1 * q[5] - 180;
	int i;
	for (i = 0; i < 6; i++)	q_DH[i] *= M_PI / 180;

	evalDHmatrix(A01, -0.5 * M_PI, .15, q_DH[0], .45, .0);
	evalDHmatrix(A12, .0, .59, q_DH[1], .0, 0);	//-.5 * pi come offset teorico ma in pratica non serve
	evalDHmatrix(A23, -.5 * M_PI, .13, q_DH[2], .0, .0);
	evalDHmatrix(A34, .5 * M_PI, .0, q_DH[3], .6471, .0);
	evalDHmatrix(A45, -.5 * M_PI, .0, q_DH[4], .0, .0);
	evalDHmatrix(A56, .0, .0, q_DH[5], .095 + offset_tool, .0);

	pose_pose_mult(A01, A12, A02);
	pose_pose_mult(A02, A23, A03);
	pose_pose_mult(A03, A34, A04);
	pose_pose_mult(A04, A45, A05);
	pose_pose_mult(A05, A56, A06);

	orientation_matrix R;

	pose_mat_getRT(A01, R, joints[0]);
	pose_mat_getRT(A02, R, joints[1]);
	pose_mat_getRT(A03, R, joints[2]);
	pose_mat_getRT(A04, R, joints[3]);
	pose_mat_getRT(A05, R, joints[4]);
	pose_mat_getRT(A06, R, joints[5]);

	//print_jointvect(q);
	/*print_posemat(A01);
	print_posemat(A02);
	print_posemat(A03);
	print_posemat(A04);
	print_posemat(A05);
	print_posemat(A06);*/
}
//calcola s, r_s dati due giunti p_a e p_b
double obst_link_distance(const translation_vector p_obs, const translation_vector p_a, const translation_vector p_b) {
	double s, dist;

	translation_vector p_ba, p_aobs, p_s, p_sum_ab, p_s_ba, p_obs_s;

	vector_diff(p_b, p_a, p_ba);
	vector_diff(p_a, p_obs, p_aobs);

	//ascissa curvilinea
	s = -vect_vect_mult(p_ba, p_aobs) / vect_vect_mult(p_ba, p_ba);
	if (s < 0)	s = 0;
	else if (s > 0)	s = 1;

	//punto corrispondente lungo il link
	scalar_vect_mult(s, p_ba, p_s_ba);
	vector_sum(p_a, p_s_ba, p_s);

	//distanza tra p_obs e p_s
	vector_diff(p_obs, p_s, p_obs_s);

	//printf("s = %6.2f\n", s);
	//print_vector(p_s);
	return vector_norm(p_obs_s);
}



typedef double quaternion_vec[4];
int get_quaternion(orientation_matrix matrix, quaternion_vec quaternion) {

	double t; //è il segno
	memset(&(quaternion[0]), 0, 4 * sizeof(double));
	quaternion[0] = .5*sqrt(matrix[0][0] + matrix[1][1] + matrix[2][2] + 1.0); //qw
	if (matrix[2][1] - matrix[1][2] >= 0) t = 1.0;
	else t = -1.0;
	quaternion[1] = .5*t*sqrt(matrix[0][0] - matrix[1][1] - matrix[2][2] + 1.0); //qx
	if (matrix[0][2] - matrix[2][0] >= 0) t = 1.0;
	else t = -1.0;
	quaternion[2] = .5*t*sqrt(matrix[1][1] - matrix[2][2] - matrix[0][0] + 1.0); //qy
	if (matrix[1][0] - matrix[0][1] >= 0) t = 1.0;
	else t = -1.0;
	quaternion[3] = .5*t*sqrt(matrix[2][2] - matrix[0][0] - matrix[1][1] + 1.0); //qz
	return 0;

}

int get_rotation(quaternion_vec quaternion, orientation_matrix matrix) {

	double t;
	memset(&(matrix[0][0]), 0, 9 * sizeof(double));
	matrix[0][0] = 2 * (quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1]) - 1;
	matrix[0][1] = 2 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3]);
	matrix[0][2] = 2 * (quaternion[1] * quaternion[3] + quaternion[0] * quaternion[2]);
	matrix[1][0] = 2 * (quaternion[1] * quaternion[2] + quaternion[0] * quaternion[3]);
	matrix[1][1] = 2 * (quaternion[0] * quaternion[0] + quaternion[2] * quaternion[2]) - 1;
	matrix[1][2] = 2 * (quaternion[2] * quaternion[3] - quaternion[0] * quaternion[1]);
	matrix[2][0] = 2 * (quaternion[1] * quaternion[3] - quaternion[0] * quaternion[2]);
	matrix[2][1] = 2 * (quaternion[2] * quaternion[3] + quaternion[0] * quaternion[1]);
	matrix[2][2] = 2 * (quaternion[0] * quaternion[0] + quaternion[3] * quaternion[3]) - 1;
	return 0;

}


// Signal handlers
void sigint_handler(int signum)
{
    stop_mod5 = true;
}
void sigabn_handler(int signum)
{
    abnormal_behaviour = true;
}


/*void error(const char *msg)
{
    perror(msg);
    //exit(0);
}*/
struct packet_in{
	double position[16];
	}__attribute__((packed));
struct packet_in cart_in; //struct centro terna utensile e orientamento

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
double mis[16];




void* Comunicazione() {
	// Timing
	struct timespec delay;
	delay.tv_sec = 0;
	delay.tv_nsec = 2E6; //2ms

	int sockfd, nbytes;
	struct sockaddr_in serv_addr;
	struct hostent *server;
	int portno = 2000;
	const char *nome_server = "192.168.0.68";
	float posiz[3];
	memset((char*)&cart_in, '0', sizeof(cart_in));


	// Socket stuff (client)
	sockfd = socket(AF_INET, SOCK_STREAM, 0);
	if (sockfd < 0) printf("ERROR opening socket");
	server = gethostbyname(nome_server);
	if (server == NULL) { fprintf(stderr, "ERROR, no such host\n"); exit(0); }
	bzero((char *)&serv_addr, sizeof(serv_addr));
	serv_addr.sin_family = AF_INET;
	bcopy((char *)server->h_addr, (char *)&serv_addr.sin_addr.s_addr, server->h_length);
	serv_addr.sin_port = htons(portno);
	if (connect(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) printf("ERROR connecting");

	float py = 3.14;

	while (1) {

		// Leggi dati dal server   
		nbytes = recv(sockfd, (void*)&cart_in, sizeof(cart_in), 0);
		if (nbytes < 0) error("ERROR reading from socket");
		//printf("centro terna:\n %f  %f  %f  \n\n", cart_in.position[0], cart_in.position[1], cart_in.position[2]);
		mis[0] = cart_in.position[0];
		mis[1] = cart_in.position[1];
		mis[2] = cart_in.position[2];
		mis[3] = cart_in.position[3];
		mis[4] = cart_in.position[4];
		mis[5] = cart_in.position[5];
		mis[6] = cart_in.position[6];
		mis[7] = cart_in.position[7];
		mis[8] = cart_in.position[8];
		mis[9] = cart_in.position[9];
		mis[10] = cart_in.position[10];
		mis[11] = cart_in.position[11];
		mis[12] = cart_in.position[12];
		mis[13] = cart_in.position[13];
		mis[14] = cart_in.position[14];
		mis[15] = cart_in.position[15];
		// Wait for the next cycle, 2 ms
		nanosleep(&delay, 0);
	}

	//chiusura socket client
	close(sockfd);
}

int main(int argc, char *argv[])
{
	const double dt = 0.002;
	int k; int i;
	int err;
	int state = 0;
	pthread_t thread;
	err = pthread_create(&thread, NULL, &Comunicazione, NULL);
	if (err != 0)
		printf("\ncan't create thread :[%s]", strerror(err));
	else
		printf("\n Thread created successfully\n");

    // Offset tool
    double z_tool = 0.19;//0.021; //non c'è più la pallina
    double x_tool = 0.0; 
   
    // Posa end effector nello spazio operativo
    double xyz_end[3];
	double phi,theta,psi; //angoli di eulero 
	double xyz_end_sim[3]; //posa end effector nello spazio operativo simulata
	//matrice di orientamento dell'end effector rispetto alla base
	double Reb[3][3];
	double Rbe[3][3]; //trasposta
	double centro_areografo_ee[3];
	double lefthandheight; //altezza mano sinistra
    // Veolocità nello spazio operativo
    double linVel[3];
    double angVel[3];
	cart_linearVel transDisp = { 0.0,0.0,0.0 };   // displacements
	cart_angularVel rotDisp = { 0.0,0.0,0.0 };    // displacements
    double t_lift;
    int enable_lift = 0;
    // posizione palline
	double KinectData[16]; //copia dati provenienti dal kinect

    int byteReceived = 0;
    bool stop_execution = false;
    bool initPacketReceived = false;
    unsigned long long loop_cnt = 0;
    int axOpenIdx = 0;
	int contatore = 0; //ausilio a kappa
	double kappa = 0.0005; //per far andare a regime la velocità dopo 4 secondi
	
	quaternion_vec errorQuat; 		//quaternione errore di orientamento EH
	orientation_matrix humanToolPoseR_T, errorR; //trasposta R_H e matrice errore orientamento EH
	
	//Aggiunti da Ale
	double temp_norm_rot; //norma omega [rad/s]
	double temp_norm;	  //norma v [m/s]
	orientation_matrix Rx_pi = {{1,0,0},{0,-1,0},{0,0,-1}}; //elementary rotation matrix of pi around x axis
	quaternion_vec errorQuat_EO; 		//quaternione errore di orientamento EO
	orientation_matrix RO, RO_rotatedT, RO_rotated, errorR_EO, errorR_EO_T; //matrice orientamento di O, trasposta R_O_rotated,  matrice errore orientamento EO e sua trasposta(per test)
	translation_vector rO;
	orientation_matrix RE_rotated, RO_T;
	orientation_matrix errorR_T;
	bool rot_sat = false; 				// flag saturazione velocità angolare 
	bool trans_sat = false;				// flag saturazione velocità lineare
	const double T_softControl = 4; 	// tempo di applicazione soft control [s];
	double k_soft = 0; 					// variabile soft control
	bool soft_start = false;			// flag che indica che va fatto il soft control
	const double tresh_position = 0.1;	// thresold for the norm of the error to shut down the control [m]
	const double tresh_rotation = 0.; 	// threshold for norm of orientation error [rad] --> 10°
	bool tran_paused = false, rot_paused = false;
	double bump_length = 0.88;	//sono tutte semi-distanze! (hp: paraurti montato centralmente)
	double bump_height = 0.25;
	double bump_depth = 0.15;
	translation_vector bump_dxH_O = {-bump_height, -bump_length, -bump_depth}; //relative a terna O
	translation_vector bump_dxL_O = {bump_height, -bump_length, -bump_depth};
	translation_vector bump_sxH_O = {-bump_height, bump_length, -bump_depth};
	translation_vector bump_sxL_O = {bump_height, bump_length, -bump_depth};
	translation_vector bump_cenH_O = {-bump_height, 0, -bump_depth};
	translation_vector bump_cenL_O = {bump_height, 0, -bump_depth};
	translation_vector bump_dxH, bump_dxL, bump_sxH, bump_sxL, bump_cenH, bump_cenL;
	double safety_cil_radius = 0.1; // raggio cilindro [m] che avvolge i link del robot per sicurezza
	double delta_safety = 0.10; // raggio cilindro [m] che avvolge i link del robot per sicurezza
	translation_vector joints[6]; //posizioni dei giunti nello spazio
	double dist[2]; //distanza ostacolo da Link1 e Link2
	double lefthand_tresh = 1.5; // threshold in [m]
	bool possible_collision = false; // flag per probabile collisione con link robot
	bool possible_collision_floor = false; //collisione paraurti suolo
	bool x_limit = false, y_limit = false, z_limit = false; //flag limiti massimi workspace raggiunti
	bool lefthand_raised = false; //flag per segnalare mano alzata sopra threshold
	bool motion_OFF = false; //disabilita ogni movimento
	translation_vector bump_points[6], bump_points_O[6]; //array di vettori posizione del paraurti
	bump_points_O = { bump_dxH_O, bump_dxL_O, bump_sxH_O, bump_sxL_O, bump_cenH_O,bump_cenL_O };
	double distances[6][2]; //matrice di distanze di ogni punto del bumper dai link 1 e 2

    // Rx and Tx packets initialization
    MVSW_DATA_PC  sx_C4GOpen_Packet_Rx;
    MVSW_DATA_PC  sx_C4GOpen_Packet_Tx;
    MVSW_SYNC_PCK sx_sync;
    MVSW_TCP_DATA UDPconn_info;

    int prev_modality = 100;
    bool isOpenControllerActive = false;
    int mod5_cnt;

    
    double t, t0;
	

    cycleTimeInfo cycletime_info;
    RTIME t1,t1Old,t2;

	joint_vector qComauJoint_ref, dqComauJoint_ref, qComauJoint_act, dqComauJoint_act;
    joint_vector qComauAxis_ref, dqComauAxis_ref, qComauAxis_act, dqComauAxis_act;
    joint_vector motorCurrent_ref;
    joint_vector qComauJoint_motion = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Create signal handlers
    signal(SIGINT, sigint_handler);
    signal(SIGABN, sigabn_handler);
    stop_mod5 = abnormal_behaviour = false;


    mlockall(MCL_CURRENT | MCL_FUTURE);


    ////////////////////////////////////////////////////////////////
    // Socket server creation and configuration                   //
    ////////////////////////////////////////////////////////////////
    UDPconn_info.si_port = C4GOPEN_SRVPORT;
    if (CreateSocketServerUDP(&UDPconn_info) == RET_ERROR)
    {
        printf("ERROR: failed to create the server socket!\n");
        exit(1);
    }
    else
    {
        UDPconn_info.si_timeout = C4GOPEN_SRVTIMEOUT;
        UDPconn_info.si_num_byte = MVZ_SYNC_PCK;
        UDPconn_info.pc_pck = NULL;
        if (!(UDPconn_info.pc_pck = (char *) malloc(UDPconn_info.si_num_byte*sizeof(char))))
        {
            printf("ERROR: failed to allocate memory for message buffer!\n");

            CloseSocketServerUDP(&UDPconn_info);
            exit(1);
        }
    }
    printf("Socket created, waiting for connection... \n");
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // RealTime task initialization                               //
    ////////////////////////////////////////////////////////////////
    RT_TASK *task = NULL;
    if (!(task = rt_task_init_schmod(nam2num(C4GCONTROL_TSK), 1, 0, 0, SCHED_FIFO, 0xFF)))
    {
        printf("ERROR: failed to create RTAI task!\n");

        CloseSocketServerUDP(&UDPconn_info);
        exit(1);
    }

    rt_task_use_fpu(task, 1);
    rt_set_oneshot_mode();
    start_rt_timer(0);
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Datalogger process creation                                //
    ////////////////////////////////////////////////////////////////
    if (!fork())
    {
        char datalog_prio[5];
        sprintf(datalog_prio, "%d", DATALOGLF_PRIO);

        execl("./datalogger", "./datalogger", argv[0], DATALOGLF_TSK, datalog_prio, DATALOGLF_MBX, NULL);
        exit(0);
    }
    printf("DATALOGGER task created.\n");

    // Wait the datalogger initialises
    while (!rt_get_adr(nam2num(DATALOGLF_TSK)))
        rt_sleep(2E6);

    // Create a mailbox
    MBX *mbx = NULL;
    if (!(mbx = rt_mbx_init(nam2num(DATALOGLF_MBX), DATALOG_MBXBUF*sizeof(char))))
    {
        printf("Cannot create mailbox to the DATALOGGER process.\n");
        exit(1);
    }
    printf("C4GCONTROL task create a mailbox to DATALOGGER.\n");
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Cycle time variable initialization                         //
    ////////////////////////////////////////////////////////////////
    cycletime_info.worstCycleTime   = 0.0;
    cycletime_info.minInterStepTime = 1000000000.0;
    cycletime_info.maxInterStepTime = -1.0;

    t1 = rt_get_time_ns();
    ////////////////////////////////////////////////////////////////

    rt_make_hard_real_time();

	// defining variables
	quaternion_vec quaternionToolPoseR, quaternionErgoPoseR; //quaternione orientamento terna utensile e quaternione orientamento terna ergonomica
	orientation_matrix humanToolPoseR, humanErgoPoseR; //matrice di orientamento terna utensile e matrice di orientamento terna ergonomica
	translation_vector humanToolPoseT, humanErgoPoseT; //posizione origine terna utensile e posizione origine terna ergonomica
	int valid = 0; //validità misure fornite dal kinect
	
	// initial values to variables
	memset(&(humanToolPoseR[0][0]), 0, 9 * sizeof(double));
	memset(&(humanToolPoseT[0]), 0, 3 * sizeof(double));
	memset(&(humanErgoPoseR[0][0]), 0, 9 * sizeof(double));
	memset(&(humanErgoPoseT[0]), 0, 3 * sizeof(double));
    // Da qui comincia la parte in tempo reale ********************************************************************************************
    while (!stop_execution && !abnormal_behaviour)
    {

	 t1Old = t1;

        // Set msg size for init packet
        if (!initPacketReceived)
            UDPconn_info.si_num_byte = MVZ_SYNC_PCK;

        // Receive and process a packet
        if ((byteReceived = ReceiveFromSocketUDP(&UDPconn_info)) > 0)
        {
            t1 = rt_get_time_ns();

            // Extract the Rx packet
            memset(&sx_C4GOpen_Packet_Rx, 0x00, sizeof(sx_C4GOpen_Packet_Rx));
            sx_C4GOpen_Packet_Rx = *(MVSW_DATA_PC*)&UDPconn_info.pc_pck[0];

            // Prepare the Tx packet
            memset(&sx_C4GOpen_Packet_Tx, 0x00, sizeof(sx_C4GOpen_Packet_Tx));
            sx_C4GOpen_Packet_Tx.NUM = sx_C4GOpen_Packet_Rx.NUM+1;
            if (sx_C4GOpen_Packet_Tx.NUM > 0x7FFF)  sx_C4GOpen_Packet_Tx.NUM = 1;

            // If a C4G_EXIT is received, stop the while loop
            if (sx_C4GOpen_Packet_Rx.AX[0].SM == C4G_EXIT)  stop_execution = true;

            // Receive the init packet and print init info
            if ((sx_C4GOpen_Packet_Rx.NUM == 10) && !initPacketReceived)
            {
                sx_sync = *(MVSW_SYNC_PCK*)&UDPconn_info.pc_pck[0];
                PrintInitializationMessage(sx_sync);

                cycletime_info.interPacketsPeriod = ((double)(sx_sync.period)) / 1000.0;

                printf("PC connected to SMP+ \n");

                initPacketReceived = true;
            }
 

            // If the robot is in DRIVE ON execute one of the open control modality
            if (sx_C4GOpen_Packet_Rx.INFPAR == C4G_STS_MOVING)
            {
                switch (sx_C4GOpen_Packet_Rx.AX[0].SM)
                {

                case C4G_DRIVING_ON:
                    // Changing modality and going into DRIVE ON...
                    if (prev_modality != C4G_DRIVING_ON)
                    {
                        prev_modality = C4G_DRIVING_ON;
                        printf("\tSoft drive on is running ... \n");
                    }

                    // Consistency check
                    if (isOpenControllerActive)  printf("INCONSISTENCY: open control active during modality 0!\n");

                    // Prepare the Tx packet
                    for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx ++)
                    {
                        sx_C4GOpen_Packet_Tx.AX[axOpenIdx].SM = C4G_DRIVING_ON;
                        sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D1 = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D1;
                        sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D2 = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D2;
                    }
                    break;


                case C4G_MOD_0:
                    // Changing modality and going into MODALITY 0...
                    if (prev_modality != C4G_MOD_0)
                    {
                        prev_modality = C4G_MOD_0;
                        printf("\tModality 0 (CLOSE) is running... \n");
                    }

                    // Consistency check
                    if (isOpenControllerActive)  printf("INCONSISTENCY: open control active during drive on!\n");

                    // Prepare the Tx packet
                    for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx ++)
                    {
                        sx_C4GOpen_Packet_Tx.AX[axOpenIdx].SM = C4G_MOD_0;
                        sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D1 = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D1;
                        sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D2 = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D2;
                    }
                    break;

                case C4G_MOD_5:
                    // Changing modality and going into MODALITY 5...
                    if (prev_modality != C4G_MOD_5)
                    {
                        prev_modality = C4G_MOD_5;
                        printf("\tModality 5 (OPEN target) is running... \n");
                    }

                    // Entering MODALITY 5...
                    if (!isOpenControllerActive)
                    {
                        isOpenControllerActive = true;

                        // Reset max position error to the actual following error
                        for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx++)
                            af_PosErrMax[axOpenIdx] = sx_sync.arm1FollowError[axOpenIdx];

                        // Reset modality 5 execution counter
                        mod5_cnt = 0;

                        // Initialise joint position for trajectory generation
                        joint_vector qComauAxis;
                        for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx ++)
                            qComauAxis[axOpenIdx] = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D1;
                        axis_to_joint_pos(qComauJoint_motion, qComauAxis);

                        // Initialise the time variables
                        t = 0.0;
                        t0 = (double)rt_get_time_ns()/1.0e9;

                    }

                    // Executing MODALITY 5...
                    if (isOpenControllerActive)
                    {
                        ///////////////////////////////////////////////////////////////////////////////////
                        // Update variables
                        ///////////////////////////////////////////////////////////////////////////////////

                        t = (double)rt_get_time_ns()/1.0e9-t0;

                        // Retrieve joint/axis/current data from C4Gopen
                        for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx ++)
                        {
                            qComauAxis_ref[axOpenIdx]   = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D1;
                            dqComauAxis_ref[axOpenIdx]  = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D2;
                            qComauAxis_act[axOpenIdx]   = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D3;
                            dqComauAxis_act[axOpenIdx]  = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D4;
                            motorCurrent_ref[axOpenIdx] = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D5;
                        }
                        axis_to_joint_pos( qComauJoint_ref, qComauAxis_ref);
                        axis_to_joint_vel(dqComauJoint_ref, dqComauAxis_ref);
                        axis_to_joint_pos( qComauJoint_act, qComauAxis_act);
                        axis_to_joint_vel(dqComauJoint_act, dqComauAxis_act);

                        // Compute the tool frame orientation wrt the base frame
                        pose_matrix tool_mpose;
                        pose_vector tool_vpose;
                        forward_kinematics_tool(qComauJoint_act, x_tool, z_tool,tool_mpose);
                        pose_mat2vect_eulZYZ(&tool_vpose, tool_mpose);

						// calcolo la tool frame sulla base del movimento simulato qComauJoint_motion
						pose_matrix tool_mpose_sim;
						pose_vector tool_vpose_sim;
						forward_kinematics_tool(qComauJoint_motion, x_tool, z_tool, tool_mpose_sim);
						pose_mat2vect_eulZYZ(&tool_vpose_sim, tool_mpose_sim);
						
						// leggo i dati provenienti dal kinect
						pthread_mutex_lock(&mutex);
						KinectData[0] = mis[0]; //x
						KinectData[1] = mis[1]; //y
						KinectData[2] = mis[2]; //z
						KinectData[3] = mis[3]; //qx
						KinectData[4] = mis[4]; //qy
						KinectData[5] = mis[5]; //qz
						KinectData[6] = mis[6]; //qw
						KinectData[7] = mis[7]; //valid
						KinectData[8] = mis[8]; //x ergo
						KinectData[9] = mis[9]; //y ergo
						KinectData[10] = mis[10]; //z ergo
						KinectData[11] = mis[11]; //qx ergo
						KinectData[12] = mis[12]; //qy ergo
						KinectData[13] = mis[13]; //qz ergo
						KinectData[14] = mis[14]; //qw ergo
						KinectData[15] = mis[15]; //left hand height
						pthread_mutex_unlock(&mutex);

						//dai dati del kinect ricavo posizione e orientamento della terna utensile (e validità di tali misure)
						humanToolPoseT[0] = KinectData[0];
						humanToolPoseT[1] = KinectData[1];
						humanToolPoseT[2] = KinectData[2];
						quaternionToolPoseR[0] = KinectData[6]; //qw
						quaternionToolPoseR[1] = KinectData[3]; //qx
						quaternionToolPoseR[2] = KinectData[4]; //qy
						quaternionToolPoseR[3] = KinectData[5]; //qz
						get_rotation(&quaternionToolPoseR[0], humanToolPoseR);
						valid = (int)KinectData[7];
						//dai dati del kinect ricavo posizione e orientamento della terna ergonomica
						humanErgoPoseT[0] = KinectData[8];
						humanErgoPoseT[1] = KinectData[9];
						humanErgoPoseT[2] = KinectData[10];
						quaternionErgoPoseR[0] = KinectData[14]; //qwergp
						quaternionErgoPoseR[1] = KinectData[11]; //qxergo
						quaternionErgoPoseR[2] = KinectData[12]; //qyergo
						quaternionErgoPoseR[3] = KinectData[13]; //qzergo
						get_rotation(&quaternionErgoPoseR[0], humanErgoPoseR);
						lefthandheight = KinectData[15]; // altezza mano sinistra


                        ///////////////////////////////////////////////////////////////////////////////////

                        ///////////////////////////////////////////////////////////////////////////////////
                        // Compute the robot motion - State Machine
                        ///////////////////////////////////////////////////////////////////////////////////
			
						joint_vector dqComauJoint = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

                        // Check that nobody during the computation raised a SIGABN
                        joint_vector dqComauAxis;
                        if (abnormal_behaviour)
                        {
                            printf("An abnormal behaviour of the system has been detected: no command to C4G allowed.\n");
                            memset(dqComauAxis, 0, 6*sizeof(double));
                        }
						else
						{

							//Start state machine ********************************************************************************************
							if (state == 0) { 	// steering robot to ergonomic pose


								//porto la posizione di O a coincidere con E (a meno di una tolleranza)
								transDisp[0] = 20.0E-2 * (humanErgoPoseT[0] - tool_vpose.pos_x) * dt;
								transDisp[1] = 20.0E-2 * (humanErgoPoseT[1] - tool_vpose.pos_y) * dt;
								transDisp[2] = 20.0E-2 * (humanErgoPoseT[2] - tool_vpose.pos_z) * dt;

								//Porto l'orientamento di O a essere uguale a quello di E (a meno di una rotazione di pi rad)
								pose_mat_getRT(tool_mpose, RO, rO); //da matrice omogenea ad mat. rotazione e traslazione

								// Errore di orientamento tra O e la terna speculare ad E
								matrix_matrix_mult(humanErgoPoseR, Rx_pi, RE_rotated); // RE_rotated = Rx_pi * RE (Rotazione di 180° della terna E)
								matrix_transpose(RO, RO_T);									// è la trasposta di RO
								matrix_matrix_mult(RE_rotated, RO_T, errorR_EO); 		// errorR = Rerg_rot * RO_T
								get_quaternion(errorR_EO, &errorQuat_EO[0]); 						// errorQuat_EO = rotm2quat(errorR_EO)

								for (i = 0; i < 3; i++) rotDisp[i] = __CONTROL_GAIN_ROT_DISPLACEMENT__ * errorQuat_EO[i + 1] * dt; // Wrob = KRw*errorQuat_EO(2:end)';*/

								//controllo su vicinananza end effector - posa ergonomica
								if (fabs(humanErgoPoseT[0] - tool_vpose.pos_x) < 0.05) transDisp[0] = 0.0;
								if (fabs(humanErgoPoseT[1] - tool_vpose.pos_y) < 0.05) transDisp[1] = 0.0;
								if (fabs(humanErgoPoseT[2] - (tool_vpose.pos_z + z_tool)) < 0.05) transDisp[2] = 0.0;

								if ((fabs(humanErgoPoseT[0] - tool_vpose.pos_x) < 0.05) &&
									(fabs(humanErgoPoseT[1] - tool_vpose.pos_y) < 0.05) &&
									(fabs(humanErgoPoseT[2] - (tool_vpose.pos_z + z_tool)) < 0.05) &&
									fabs(errorQuat_EO[1]) < 0.05 && fabs(errorQuat_EO[2]) < 0.05 && fabs(errorQuat_EO[3]) < 0.05) 
								{
									state = 1;
									printf("**********************\n**********************\nGO TO STATE 1!\n**********************\n**********************");
								}
								
							}else if (state == 1) { 	// wait for human gesture to start

								transDisp[0] = transDisp[1] = transDisp[2] = 0.0;
								rotDisp[0] = rotDisp[1] = rotDisp[2] = 0.0;

								if (valid == 1 && lefthandheight > 1.5) { //se alzo la mano sinistra quando il paraurti è nella posizione ergonomica allora attivo il controllo
									state = 2;
									printf("**********************\n**********************\nGO TO STATE 2!\n**********************\n**********************");
								}

							}
							else if (state == 2) {		// computing control law (robot ee velocities)

								double t_al = 1.5; //tempo di accelerazione
								double vel_l, acc_l;
								double t_inizl = t_lift, Tl = 15;
								double dir_vel_l[3] = { 0.0 ,0.0, 1 };

								vel_l = 0.005;
								acc_l = vel_l / t_al;
								double freq = 0.10; //frequenza velocità profilo sinusoidale

								i = 0;

								//soft start
								/*if (contatore < 2000) kappa = 0.0005*contatore;
								else kappa = 1;*/

								//controllo traslazione
								//for (i = 0; i < 3; i++) transDisp[i] = __CONTROL_GAIN_TRANS_DISPLACEMENT__ * (humanErgoPoseT[i] - humanToolPoseT[i])*dt; // Vrob = KRv*(Perg - Phum);
								transDisp[0] = kappa * __CONTROL_GAIN_TRANS_DISPLACEMENT__ * (humanErgoPoseT[0] - humanToolPoseT[0]) * dt;
								transDisp[1] = kappa * __CONTROL_GAIN_TRANS_DISPLACEMENT__ * (humanErgoPoseT[1] - humanToolPoseT[1]) * dt;
								transDisp[2] = kappa * __CONTROL_GAIN_TRANS_DISPLACEMENT__ * (humanErgoPoseT[2] - humanToolPoseT[2]) * dt;

								//Controllo rotazionale
								matrix_transpose(humanToolPoseR, humanToolPoseR_T);// è la trasposta della matrice di orientamento del tool frame rispetto al world frame
								matrix_matrix_mult(humanErgoPoseR, humanToolPoseR_T, errorR); // errorR = Rerg*Rhum';
								get_quaternion(errorR, &errorQuat[0]); // errorQuat = rotm2quat(errorR);

								//test errore di orientamento trasposto!
								//matrix_transpose(&errorR[0][0], &errorR_T[0][0]);
								//get_quaternion(&errorR_T[0][0], &errorQuat[0]);

								for (i = 0; i < 3; i++) rotDisp[i] = __CONTROL_GAIN_ROT_DISPLACEMENT__ * errorQuat[i + 1] * dt; // Wrob = KRw*errorQuat(2:end)';*/

								//update contatore mod 5
								contatore = contatore + 1;


								/*// traiettoria sinusoidale
								if (t<10)
								{
									transDisp[0] = transDisp[1] = transDisp[2] = 0.0;
									rotDisp[0] = rotDisp[1] = rotDisp[2] = 0.0;
								}
								else if ((t >= 10) && (t < 50))
								{
									transDisp[0] = 0.0;
									transDisp[1] = 0.0;
									transDisp[2] = 0.1*cos(2 * 3.14*freq*t)*dt;
									rotDisp[0] = rotDisp[1] = rotDisp[2] = 0.0;
								}
								else {
									transDisp[0] = transDisp[1] = transDisp[2] = 0;
									rotDisp[0] = rotDisp[1] = rotDisp[2] = 0;
								}*/


							}
							else { // idle state!
								transDisp[0] = 0.0;
								transDisp[1] = 0.0;
								transDisp[2] = 0.0;
								rotDisp[0] = 0.0;
								rotDisp[1] = 0.0;
								rotDisp[2] = 0.0;
							}
							//End state machine ********************************************************************************************

							///////////////////////////////////////////////////////////////////////////////////
							// Check the computed motion - SAFETY
							///////////////////////////////////////////////////////////////////////////////////
							//reset flags
							motion_OFF = false;
							possible_collision = false;
							possible_collision_floor = false;
							lefthand_raised = false;
							x_limit = false;
							y_limit = false;
							z_limit = false;
							trans_sat = false;
							rot_sat = false;
							
							// Left hand height
							if (lefthandheight > lefthand_tresh)	lefthand_raised = true;

							//Limiti end effector in terna W
							if ((tool_vpose.pos_x < 0.55) && (transDisp[0] < 0) || (tool_vpose.pos_x > 1.00) && (transDisp[0] > 0)) {
								transDisp[0] = 0.0;
								x_limit = true;
							}
							if ((tool_vpose.pos_y < -1.00) && (transDisp[1] < 0) || (tool_vpose.pos_y > 1.00) && (transDisp[1] > 0)) {
								transDisp[1] = 0.0;
								y_limit = true;
							}
							if ((tool_vpose.pos_z+z_tool < bump_length) && (transDisp[2] < 0) || (tool_vpose.pos_z+z_tool > 1.30) && (transDisp[2] > 0)) {
								transDisp[2] = 0.0;
								z_limit = true;
							}

							//saturazione velocità massima traslazione
							temp_norm = vector_norm(transDisp);
							if (temp_norm > __TRAN_MAX__) {
								for (i = 0; i < 3; i++) transDisp[i] = transDisp[i] * (__TRAN_MAX__ / temp_norm);
								trans_sat = true;
							}

							//saturazione massima velocità angolare
							temp_norm_rot = vector_norm(rotDisp);
							if (temp_norm_rot > __ROT_MAX__) {
								for (i = 0; i < 3; i++) rotDisp[i] = rotDisp[i] * (__ROT_MAX__ / temp_norm_rot);
								rot_sat = true;
							}

							// alternative computation of forward kinematics
							fkine(qComauJoint_act, joints);

							//Controllo collisioni paraurti contro robot o suolo
							/*pose_vector_mult(tool_mpose, bump_dxH_O, bump_dxH);
							pose_vector_mult(tool_mpose, bump_dxL_O, bump_dxL);
							pose_vector_mult(tool_mpose, bump_sxH_O, bump_sxH);
							pose_vector_mult(tool_mpose, bump_sxL_O, bump_sxL);
							pose_vector_mult(tool_mpose, bump_cenH_O, bump_cenH);
							pose_vector_mult(tool_mpose, bump_cenL_O, bump_cenL);*/
							for (i = 0;i < 6 && !possible_collision && !possible_collision_floor; i++) {
								pose_vector_mult(tool_mpose, bump_points_O[i], bump_points[i]);
								dist[i][0] = obst_link_distance(bump_points[i], joints[0], joints[1]);
								dist[i][1] = obst_link_distance(bump_points[i], joints[2], joints[3]);

								//confronto con soglie sicurezza link
								if (dist[i][0] < delta_safety || dist[i][1] < delta_safety) possible_collision = true;

								//confronto altezza da suolo
								if (bump_points[i][2] < 0.1)	possible_collision_floor = true; //10 [cm] di soglia
							}

							//calcolo distanze tra ogni punto del bumper e i link comau (Link1 = giunti 0-1; Link2 = giunti 2-3)
							/*dist[0] = obst_link_distance(bump_dxH, joints[0], joints[1]);
							dist[1] = obst_link_distance(bump_dxH, joints[2], joints[3]);*/
					

							/*
							// check dead-zones
							if (fabs(humanErgoPoseT[0] - humanToolPoseT[0]) < tresh_position &&
								fabs(humanErgoPoseT[1] - humanToolPoseT[1]) < tresh_position &&
								fabs(humanErgoPoseT[2] - humanToolPoseT[2]) < tresh_position){
								tran_paused = true;
							}else{
								tran_paused = false;
								//contatore = 0;
							}
							if(fabs(errorQuat[1]) < tresh_rotation &&
								fabs(errorQuat[2]) < tresh_rotation &&
								fabs(errorQuat[3]) < tresh_rotation){
								rot_paused = true;
							}else{
								rot_paused = false;
								//contatore = 0;
							}
							*/

							// stops every movement of the robot
							if (valid != 1 || lefthand_raised || possible_collision) {
								transDisp[0] = 0.0;
								transDisp[1] = 0.0;
								transDisp[2] = 0.0;
								rotDisp[0] = 0.0;
								rotDisp[1] = 0.0;
								rotDisp[2] = 0.0;
							}
						}

						///////////////////////////////////////////////////////////////////////////////////
						// Calculate the joint displacement - Prepare Packet - Print e Log info
						///////////////////////////////////////////////////////////////////////////////////

                        joint_vector dJoint_motion;
                        if (inverse_diffkinematics_tool(qComauJoint_act, dJoint_motion, transDisp, rotDisp, 0.0, z_tool)<0){
                            printf("Singularity in diffkinematics inversion.\n");
                            abnormal_behaviour = true;
                        }
                        else
                        {
                            int k;
                            for (k=0; k<6; k++)
                                 qComauJoint_motion[k] += dJoint_motion[k];
                        }
                        // solo per la simulazione --> in questo modo riesco a vedere la simulazione, ma posso anche applicare i comandi al comau
						/* if (inverse_diffkinematics_tool(qComauJoint_motion, dJoint_motion, transDisp, rotDisp, 0.0, z_tool) < 0) {
							printf("Singularity in diffkinematics inversion. \n");
							abnormal_behaviour = true;
						}
						else {
							int k;
							for (k = 0; k < 6; k++)
								qComauJoint_motion[k] += dJoint_motion[k];
						} */
                      
                        k = 0;

                        // Calcola il riferimento incrementale
                        for (k=0; k<6; k++)
                        {
                             dqComauJoint[k] = qComauJoint_motion[k] - qComauJoint_ref[k];
                        }
                        joint_to_axis_vel(dqComauJoint, dqComauAxis);



						//Manipolazione dati solo per print e datalog (non utili ai fini dell'algo)
						// Lettura cinematica diretta
                        xyz_end[0] = tool_vpose.pos_x;
                        xyz_end[1] = tool_vpose.pos_y;
                        xyz_end[2] = tool_vpose.pos_z; 
						//lettura cinematica diretta simulata
						xyz_end_sim[0] = tool_vpose_sim.pos_x;
						xyz_end_sim[1] = tool_vpose_sim.pos_y;
						xyz_end_sim[2] = tool_vpose_sim.pos_z;
						//Lettura angoli di eulero
						phi = tool_vpose.phi;
						theta = tool_vpose.theta;
						psi = tool_vpose.psi;

						Reb[0][0] = cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi);
						Reb[0][1] = -cos(phi)*cos(theta)*sin(psi)-sin(phi)*cos(psi);
						Reb[0][2] = cos(phi)*sin(theta);
						Reb[1][0] = sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi);
						Reb[1][1] = -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi);
						Reb[1][2] = sin(phi)*sin(theta);
						Reb[2][0] = -sin(theta)*cos(psi);
						Reb[2][1] = sin(theta)*sin(psi);
						Reb[2][2] = cos(theta);
						//calcolo posizione areografo rispetto all'end effector
						matrix_transpose(Reb, Rbe); //trasposta della rotazione
						double diff[3];
						diff[0] = humanToolPoseT[0] - xyz_end[0];
						diff[1] = humanToolPoseT[1] - xyz_end[1];
						diff[2] = humanToolPoseT[2] - xyz_end[2];
						matrix_vector_mult(Rbe, diff, centro_areografo_ee);

                        // Riferimento nello spazio dei giunti
                        double ref_cart;
                        ref_cart = linVel[2];

						// Display info on screen
                        if ((mod5_cnt % 250)==0) //print info every 0.5 s
                        {
                            printf("\n\ntime: %.3f *******************\nState = %d\n", t, state);
                            //printf("dq_cont (incremento motore calcolato):  %.3e   %.3e   %.3e   %.3e   %.3e   %.3e\n",dqComauAxis[0],dqComauAxis[1],dqComauAxis[2],dqComauAxis[3],dqComauAxis[4],dqComauAxis[5]);
							//printf("qComauJoint_ref(riferimento giunti):  %.3e   %.3e   %.3e   %.3e   %.3e   %.3e\n",qComauJoint_ref[0],qComauJoint_ref[1],qComauJoint_ref[2],qComauJoint_ref[3],qComauJoint_ref[4],qComauJoint_ref[5]);
							//printf("qComauJoint_act (attuale):      %.3e   %.3e   %.3e   %.3e   %.3e   %.3e\n\n",qComauJoint_act[0],qComauJoint_act[1],qComauJoint_act[2],qComauJoint_act[3],qComauJoint_act[4],qComauJoint_act[5]);
							//printf("qComauJoint_motion (simulato):       %.3e   %.3e   %.3e   %.3e   %.3e   %.3e\n\n", qComauJoint_motion[0], qComauJoint_motion[1], qComauJoint_motion[2], qComauJoint_motion[3], qComauJoint_motion[4], qComauJoint_motion[5]);
							//printf("tool position:  %.3e  %.3e  %.3e\n\n", humanToolPoseT[0], humanToolPoseT[1], humanToolPoseT[2]);
							//printf("ergo position:  %.3e  %.3e  %.3e\n\n", humanErgoPoseT[0], humanErgoPoseT[1], humanErgoPoseT[2]);
							//printf("dJoint_motion (simulato):  %.3e   %.3e   %.3e   %.3e   %.3e   %.3e\n\n", dJoint_motion[0], dJoint_motion[1], dJoint_motion[2], dJoint_motion[3], dJoint_motion[4], dJoint_motion[5]);
							printf("transDisp end effector:  %.3e  %.3e  %.3e\n", transDisp[0], transDisp[1], transDisp[2]);
							printf("rotDisp end effector:  %.3e  %.3e  %.3e\n", rotDisp[0], rotDisp[1], rotDisp[2]);
							//printf("END EFFECTOR Simulato:  %.3e  %.3e  %.3e\n\n", xyz_end_sim[0], xyz_end_sim[1], xyz_end_sim[2]);
							printf("END EFFECTOR VERO: %.3e  %.3e  %.3e\n", xyz_end[0], xyz_end[1], xyz_end[2]);
							//printf("aerografo rispetto ad end effector: %.3e  %.3e  %.3e\n\n", centro_areografo_ee[0], centro_areografo_ee[1], centro_areografo_ee[2]);
							printf("altezza mano sinistra:  %.3e\n", lefthandheight);
							//printf("error_Quat_EO: %.3e  %.3e  %.3e\n", errorQuat_EO[1], errorQuat_EO[2], errorQuat_EO[3]);
							printf("norma rotazione = %.3e ; rotation saturated = %d \n", temp_norm_rot, rot_sat);
							//printf("error_Quat: %.3e  %.3e  %.3e\n", errorQuat[1], errorQuat[2], errorQuat[3]);
							printf("tran_paused = %d; rot_paused = %d\n", tran_paused, rot_paused);
							printf("possible_collision = %d\n", possible_collision);
							printf("x_limit = %d;	y_limit = %d;	z_limit = %d;\n",x_limit,y_limit,z_limit);
							printf("estremoDX bumper: (%.3e  %.3e  %.3e)\n", bump_dxH[0], bump_dxH[1], bump_dxH[2]);
							printf("estremoSX bumper: (%.3e  %.3e  %.3e)\n", bump_sxH[0], bump_sxH[1], bump_sxH[2]);
							//for (i = 0;i < 6;i++)	printf("joint %i: (%6.3f, %6.3f, %6.3f)\n", i + 1, joints[i][0], joints[i][1], joints[i][2]);
							printf("Collisione con links = %d;	collisione suolo = %d\n",possible_collision, possible_collision_floor);
                        }

						// Preparazione pacchetto da inviare a C4GOpen
                        for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx++)
                        {
                            sx_C4GOpen_Packet_Tx.AX[axOpenIdx].SM = C4G_MOD_5;

                            // Update open axis positions and velocities

                            // NON SCRIVE
                           //sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D1 = qComauAxis_ref[axOpenIdx];
                           //sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D2 = 0.0;

                            // SCRIVE
                            sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D1 = qComauAxis_ref[axOpenIdx] + dqComauAxis[axOpenIdx];
                            sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D2 = dqComauAxis[axOpenIdx];

                            // Update the following error and check weather it is within the specified bounds
                            af_PosErr[axOpenIdx] = (sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D1 - sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D3);



                            if (af_PosErr[axOpenIdx] < 0)
                            {
                                if (af_PosErr[axOpenIdx] < (-af_PosErrMax[axOpenIdx]))
                                    sx_C4GOpen_Packet_Tx.AX[axOpenIdx].SM = C4G_FOLLOWING_ERROR;
                            }
                            else
                            {
                                if (af_PosErr[axOpenIdx] > ( af_PosErrMax[axOpenIdx]))
                                    sx_C4GOpen_Packet_Tx.AX[axOpenIdx].SM = C4G_FOLLOWING_ERROR;
                            }
                        }

                        // Send data to datalogger
                        datalog_msg msg;
                        msg.cmd = DATALOG_NONE;
                        msg.numVar = 115;
                        msg.varVector[0] = t;

                        int k;
                        for (k=0; k<6; k++)
                        {
                            msg.varVector[ 1+k] = qComauJoint_ref[k];
                            msg.varVector[ 7+k] = qComauJoint_act[k];
                            msg.varVector[13+k] = qComauAxis_ref[k];
                            msg.varVector[19+k] = qComauAxis_act[k];
                            msg.varVector[25+k] = dqComauJoint_ref[k];
                            msg.varVector[31+k] = dqComauJoint_act[k];
                            msg.varVector[37+k] = dqComauAxis_ref[k];
                            msg.varVector[43+k] = dqComauAxis_act[k];
							
							msg.varVector[78+k] = qComauJoint_motion[k];
							//for(i=0;i<3;i++)	msg.varVector[96+k+i] = joints[k][i];
							

                        }
                        msg.varVector[49] = tool_vpose_sim.pos_x; //posa simulata dell'end effector
                        msg.varVector[50] = tool_vpose_sim.pos_y;
                        msg.varVector[51] = tool_vpose_sim.pos_z;
						msg.varVector[52] = transDisp[0]; //incrementi di posizione dell'end effector (velocità)
                        msg.varVector[53] = transDisp[1];
                        msg.varVector[54] = transDisp[2];
						msg.varVector[55] = humanToolPoseT[0]; //inserisco i dati dal kinect nel datalogger
						msg.varVector[56] = humanToolPoseT[1];
						msg.varVector[57] = humanToolPoseT[2];
						msg.varVector[58] = valid;
						msg.varVector[59] = humanErgoPoseT[0];
						msg.varVector[60] = humanErgoPoseT[1];
						msg.varVector[61] = humanErgoPoseT[2];
						msg.varVector[62] = tool_vpose.pos_x;//posa reale end effector
						msg.varVector[63] = tool_vpose.pos_y;
						msg.varVector[64] = tool_vpose.pos_z;
						msg.varVector[65] = centro_areografo_ee[0];
						msg.varVector[66] = centro_areografo_ee[1];
						msg.varVector[67] = centro_areografo_ee[2];
						msg.varVector[68] = quaternionErgoPoseR[0];
						msg.varVector[69] = quaternionErgoPoseR[1];
						msg.varVector[70] = quaternionErgoPoseR[2];
						msg.varVector[71] = quaternionErgoPoseR[3];
						msg.varVector[72] = quaternionToolPoseR[0];
						msg.varVector[73] = quaternionToolPoseR[1];
						msg.varVector[74] = quaternionToolPoseR[2];
						msg.varVector[75] = quaternionToolPoseR[3];
						msg.varVector[76] = tran_paused;
						msg.varVector[77] = rot_paused;
						msg.varVector[84] = state;
						msg.varVector[85] = lefthandheight;
						msg.varVector[86] = errorQuat[0];
						msg.varVector[87] = errorQuat[1];
						msg.varVector[88] = errorQuat[2];
						msg.varVector[89] = errorQuat[3];
						msg.varVector[90] = __CONTROL_GAIN_TRANS_DISPLACEMENT__;
						msg.varVector[91] = __CONTROL_GAIN_ROT_DISPLACEMENT__;
						msg.varVector[92] = possible_collision;
						msg.varVector[93] = x_limit;
						msg.varVector[94] = y_limit;
						msg.varVector[95] = z_limit;
						msg.varVector[96] = lefthand_raised;
						msg.varVector[97] = possible_collision_floor;


                        if (rt_mbx_send_if(mbx, &msg, sizeof(msg)) < 0)
                            printf("Error sending message to datalogger.\n");
						///////////////////////////////////////////////////////////////////////////////////

						// Update modality 5 counter
                        mod5_cnt++;

                    }// end if(IsOpenControllerActive)

					// Exiting MODALITY 5...
                    if (stop_mod5 || abnormal_behaviour)
                    {
                        isOpenControllerActive = false;

                        // Send the EXIT FROM OPEN command
                        sx_C4GOpen_Packet_Tx.INFPAR = C4G_EXIT_FROM_OPEN;
                        for (axOpenIdx = 0; axOpenIdx < sx_sync.nAxOpen; axOpenIdx++)
                        {
                            sx_C4GOpen_Packet_Tx.AX[axOpenIdx].SM = C4G_MOD_5;
                            sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D1 = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D1;
                            sx_C4GOpen_Packet_Tx.AX[axOpenIdx].D2 = sx_C4GOpen_Packet_Rx.AX[axOpenIdx].D2;
                        }

                        printf("\t... exiting modality 5 (OPEN target)\n");
                    }
                    break;

					default:
                    printf("\tModality %d not available... \n", (int)sx_C4GOpen_Packet_Rx.AX[0].SM);
                    break;
                }
            }// end of "if the robot is in DRIVE ON"
            else    // The robot is in DRIVE OFF
            {
                // Echo the received message
                memcpy(&sx_C4GOpen_Packet_Tx,&sx_C4GOpen_Packet_Rx,byteReceived);
                sx_C4GOpen_Packet_Tx.NUM +=1;
                if (sx_C4GOpen_Packet_Tx.NUM >0x7FFF)  sx_C4GOpen_Packet_Tx.NUM = 1;
            }  // End if STS_MOVING



            ////////////////////////////////////////////////////////////////
            // Send the Tx message to the C4G                             //
            ////////////////////////////////////////////////////////////////
            UDPconn_info.si_num_byte = byteReceived;
            memcpy(UDPconn_info.pc_pck,(unsigned char*)&sx_C4GOpen_Packet_Tx, UDPconn_info.si_num_byte);

            if (SendFromSocketUDP(&UDPconn_info) < byteReceived)
                printf("ERROR: too few bytes have been sent!");
            ////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////
            // Monitor the cycle time                                     //
            ////////////////////////////////////////////////////////////////
            if (loop_cnt >= 15)
            {
                t2 = rt_get_time_ns();

                // Evaluate the time required to reply to the SMP+ message
                cycletime_info.lastCycleTime = (double)(t2 - t1) / 1000000000.0;

                // If this time exceeds the inter-packet period terminate the execution
                if (cycletime_info.lastCycleTime > cycletime_info.interPacketsPeriod)
                {
                    printf("SCHEDULING ERROR: cycle time exceeded C4G <-> PC inter-packets period! (t1:%.9f t2: %.9f)\n",(double)(t1) / 1000000000.0,(double)(t2) / 1000000000.0);
                    stop_execution = true;
                }

                // Update the cycle time statistics
                if (cycletime_info.lastCycleTime > cycletime_info.worstCycleTime) cycletime_info.worstCycleTime = cycletime_info.lastCycleTime;

                // Evaluate the inter-step time
                cycletime_info.lastInterStepTime = ((double)(t1 - t1Old)) / 1000000000.0;
                if (cycletime_info.lastInterStepTime > cycletime_info.maxInterStepTime)  cycletime_info.maxInterStepTime = cycletime_info.lastInterStepTime;
                if (cycletime_info.lastInterStepTime < cycletime_info.minInterStepTime)  cycletime_info.minInterStepTime = cycletime_info.lastInterStepTime;
            }
            ////////////////////////////////////////////////////////////////

            loop_cnt++;
        }// end of if ((byteReceived = ReceiveFromSocketUDP(&UDPconn_info)) > 0)
    } // end while


    ////////////////////////////////////////////////////////////////
    // Close the socket server                                    //
    ////////////////////////////////////////////////////////////////
    CloseSocketServerUDP(&UDPconn_info);
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Terminate the datalogger process                           //
    ////////////////////////////////////////////////////////////////
    // Sending the quit command to datalogger
    datalog_msg quit_msg;
    quit_msg.cmd = DATALOG_QUIT;
    quit_msg.numVar = 0;

    if (rt_mbx_send_if(mbx, &quit_msg, sizeof(quit_msg)) < 0)
        printf("C4GCONTROL task: error sending quit message to mbx.\n");

    // Waiting that datalogger ends
    while (rt_get_adr(nam2num(DATALOGLF_TSK)))
        rt_sleep(nano2count(1E9));

    // Delete the mailbox to datalogger
    rt_mbx_delete(mbx);
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Delete the Real Time task                                  //
    ////////////////////////////////////////////////////////////////
    rt_make_soft_real_time();
    stop_rt_timer();
    rt_task_delete(task);
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Compute and display cycle time statistics                  //
    ////////////////////////////////////////////////////////////////
    if (loop_cnt > 15)  PrintCycleTimeStatistics(cycletime_info);
    ////////////////////////////////////////////////////////////////

    return 0;
}
