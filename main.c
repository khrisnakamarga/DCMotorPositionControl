/* Lab #8 - Khrisna Kamarga */

/* includes */
#include "stdio.h"
#include "MyRio.h"
#include "me477.h"
#include <string.h>
#include <pthread.h>
#include "TimerIRQ.h" // TimerIRQ interrupt
#include "IRQConfigure.h"
#include "matlabfiles.h"
#include "math.h"
#include "ctable2.h"

// For PIDF header
struct biquad {
	double b0; double b1; double b2; // numerator
	double a0; double a1; double a2; // denominator
	double x0; double x1; double x2; // input
	double y1; double y2; }; // output

#include "myPIDF.h"

#define SATURATE(x,lo,hi) ((x) < (lo) ? (lo) : (x) > (hi) ? (hi) : (x))

typedef struct {double xfa; // position
				double v; // velocity limit
				double a; // acceleration limit
				double d; // dwell time (s)
} seg;

typedef struct {
	NiFpga_IrqContext irqContext; // IRQ context reserved
	struct table *a_table; // table
	seg *profile; // profile
	int nseg; // no. of segs
	NiFpga_Bool irqThreadRdy; // IRQ thread ready flag
} ThreadResource;

struct table {
	const char *e_label; // entry label label
	int e_type; // entry type (0-show; 1-edit
	double value; // value
};

/* prototypes */

void *Timer_Irq_Thread(void* resource);

double cascade(double xin, // input
		struct biquad *fa, // biquad array
		int ns, // no. segments
		double ymin, // min output
		double ymax); // max output

double pos(void);

int Sramps(seg *segs, int *iseg, int nseg, int *itime, double T, double *xa);

/* definitions */

MyRio_Aio CI0;
MyRio_Aio CO0;
uint32_t timeoutValue;
double VADin;
double VADout;
NiFpga_Session myrio_session;
MyRio_Encoder encC0;
NiFpga_Session myrio_session;
double e; // current error
double hisat = 7.5; // Max saturation voltage
double losat = -7.5; // Min saturation value
static double k = 0.0451; // constant k = kvi * kt (kvi = 0.41, kt = 0.11)

// MATLAB Data
#define IMAX 4000 // # of max points = 500
static double pact_buffer[IMAX]; // pact buffer (250 pts)
static double pref_buffer[IMAX]; // pref buffer (250 pts)
static double torq_buffer[IMAX]; // torque buffer (250 pts)
static double *pact_pt = pact_buffer; // pact buffer pointer
static double *torq_pt = torq_buffer; // torque buffer pointer
static double *pref_pt = pref_buffer; // pref buffer pointer
MATFILE *mf;
int err;
// MATLAB Datafct

/*---------------------------------------------------------------------
Function main()
Purpose: creates the timer interrupt thread that performs
		 real-time DC motor position control. Calls the
		 ctable2 function to display the measured values.
		 The user can press the DEL key to stop the motor.
*---------------------------------------------------------------------*/
int main(int argc, char **argv) {
	MyRio_IrqTimer irqTimer0;
	ThreadResource irqThread0;
	pthread_t thread;
	NiFpga_Status status;
	int nseg;

	double vmax = 50.; // revs/s
	double amax = 20.; // revs/s^2
	double dwell = 1.0; // s
	seg mySegs[8] = { // revs
		{10.125, vmax, amax, dwell},
		{20.250, vmax, amax, dwell},
		{30.375, vmax, amax, dwell},
		{40.500, vmax, amax, dwell},
		{30.625, vmax, amax, dwell},
		{20.750, vmax, amax, dwell},
		{10.875, vmax, amax, dwell},
		{ 0.000, vmax, amax, dwell}
	};
	nseg = 8;



	/* Initialize Table Editor Variables */
	char *Table_Title = "Motor Control Table";
	table my_table[] = {
		{"P_ref: (rev) ", 1, 0.0 },
		{"P_act: (rev) ", 0, 0.0 },
		{"VDAout: (mV) ", 0, 0.0 }
	};

    status = MyRio_Open();		    /*Open the myRIO NiFpga Session.*/
    if (MyRio_IsNotSuccess(status)) return status;

    //my code here
    // Registers corresponding to the IRQ channel
    irqTimer0.timerWrite = IRQTIMERWRITE;
    irqTimer0.timerSet = IRQTIMERSETTIME;
    timeoutValue = 5000; // BTI = 5ms

    // configure TimerIRQ
    status = Irq_RegisterTimerIrq(
    			&irqTimer0,
				&irqThread0.irqContext,
				timeoutValue);

    // Set the indicator to allow the new thread
    irqThread0.irqThreadRdy = NiFpga_True;

    /* Set table to irqthread table */
    irqThread0.a_table = my_table;
    irqThread0.profile = mySegs;
    irqThread0.nseg = nseg;

    // Create new thread to catch the IRQ
    /* Create New Thread to Catch the IRQ */
    status = pthread_create(&thread,
							NULL,
							Timer_Irq_Thread,
							&irqThread0);

    /* Run table editor */
    ctable2(Table_Title, my_table, 6);

    // Signal the new thread to terminate
    irqThread0.irqThreadRdy = NiFpga_False;
    status = pthread_join(thread, NULL);

    // Unregistering the timer interrupt
    status = Irq_UnregisterTimerIrq( &irqTimer0,
    								 irqThread0.irqContext);

	status = MyRio_Close();	 /*Close the myRIO NiFpga Session. */
	return status;
}

/*---------------------------------------------------------------------
Function Timer_Irq_Thread()
Purpose: Performs real-time DC motor velocity control by using
		 the transfer function coefficients defined in main()
Parameter:
	void* resource = pointer to the ThreadResource for the ISR
*---------------------------------------------------------------------*/
void *Timer_Irq_Thread(void* resource) {
	extern NiFpga_Session myrio_session; // getting myrio session from main()
	uint32_t irqAssert;
	// cast thread resource to ThreadResource structure
	ThreadResource* threadResource = (ThreadResource*) resource;

	/* Declare Names for Table Entries From Table Pointer */
	double *pref = &((threadResource->a_table+0)->value);
	double *pact = &((threadResource->a_table+1)->value);
	double *VDAout = &((threadResource->a_table+2)->value);
	seg *mySegs = threadResource->profile;
	int nseg = threadResource->nseg;
	int itime = -1;
	int iseg = -1;
	uint32_t timeoutValue = 5000; // timeout value in microseconds
	double T = 5.0 / 1000.0; // 5 ms;

	// initialize analog interfaces before allowing IRQ
	AIO_initialize(&CI0, &CO0); // initialize
	Aio_Write(&CO0, 0.0); // zero analog input
	EncoderC_initialize(myrio_session, &encC0);


	// initialize Sramps
	Sramps(mySegs, &iseg, nseg, &itime, T, pref);


	irqAssert = 0;

	while (threadResource->irqThreadRdy == NiFpga_True) {
		// wait for the occurrence (or timeout) of the IRQ.
		Irq_Wait( threadResource->irqContext,
				  TIMERIRQNO,
				  &irqAssert,
				  (NiFpga_Bool*) &(threadResource->irqThreadRdy));
		// schedule the next interrupt
		NiFpga_WriteU32( myrio_session,
						 IRQTIMERWRITE,
						 timeoutValue );
		NiFpga_WriteBool( myrio_session,
						  IRQTIMERSETTIME,
						  NiFpga_True);
		// if there is an occurrence of the IRQ
		if (irqAssert) {
			//
			Sramps(mySegs, &iseg, nseg, &itime, T, pref);
			*pact = pos() * (1/2000.0); // revolution
			e = (*pref - *pact)*(2*M_PI);
			*VDAout = cascade(e, PIDF, PIDF_ns, losat, hisat);

			Aio_Write(&CO0, *VDAout); // debugging (send y(n) to AO)
			*VDAout = *VDAout * 1000;

			// MATLAB Data
			if (pact_pt < pact_buffer + IMAX) {
				*torq_pt++ = (*VDAout) * k; // mNm
				*pact_pt++ = *pact;
				*pref_pt++ = *pref;
			}

			Irq_Acknowledge(irqAssert); // acknowledge the occurrence
		}
	}

	double bti = T*1000.0;
	/* Comment out the below section when debugging */
	mf = openmatfile("Lab8_Khrisna.mat", &err);
	if(!mf) printf("Can't open mat file %d\n", err);
	matfile_addstring(mf, "myName", "Khrisna Kamarga");
	matfile_addmatrix(mf, "pact", pact_buffer, IMAX, 1, 0);
	matfile_addmatrix(mf, "torq", torq_buffer, IMAX, 1, 0);
	matfile_addmatrix(mf, "pref", pref_buffer, IMAX, 1, 0);
	matfile_addmatrix(mf, "pidf", (double *) PIDF, 1, 11, 0);
	matfile_addmatrix(mf, "bti", &bti, 1, 1, 0);
	matfile_close(mf);
	/* End of MATLAB Section */

	// terminate the thread and return from function
	pthread_exit(NULL);
	return NULL;
}

/*---------------------------------------------------------------------
Function cascade()
Purpose: Calculates the signal output based on Biquad approximation
Parameter:
	double xin = input signal
	struct biquad *fa = pointer to the biquad coefficients struct
	int ns = number of biquad segments
	double ymin = minimum saturation output
	double ymax = maximum saturation output
*---------------------------------------------------------------------*/
double cascade( double xin, // input
				struct biquad *fa, // biquad array
				int ns, // no. segments
				double ymin, // min output
				double ymax) { // max output
	struct biquad* f;
	double y0;
	f = fa; // pointer to the first biquad
	y0 = xin; // input as output

	f->x0 = y0; // y0 as input value becomes x0 of the current biquad
	// biquad cascade function
	y0 = (f->b0*f->x0 + f->b1*f->x1 + f->b2*f->x2 - f->a1*f->y1 - f->a2*f->y2)
		  /f->a0;

	y0 = SATURATE(y0, ymin, ymax); // saturate y0
	// update the previous value of x and y as next
	f->x2 = f->x1;
	f->x1 = f->x0;
	f->y2 = f->y1;
	f->y1 = y0;
	f++; // increment pointer to get the next biquad
	return y0;
}

/*--------------------------------------------------------------
Function pos()
Purpose: Calculating the position of the motor by accessing the
		 encoder counter with the knowledge of the BDI
--------------------------------------------------------------*/
double pos(void) {
	static int32_t cn; // current encoder count
	static int32_t cn1; // previous encoder count
	static int count = 1; // number of times vel() is called
	double encoderPos; // cn - cn1
	cn = Encoder_Counter(&encC0); // reads the encoder counter
	if (count == 1) {
		cn1 = cn; // initial condition is set the same as the first count
		count = 0; // ensures that cn1 = cn is only done on the first call
	}
	encoderPos= cn - cn1; // current and previous count difference
	return encoderPos; // return the rpm (double)
}

/*
 * Sramps.c
 *
 *  Created on: Mar 18, 2016
 *      Author: garbini
 */

int Sramps(seg *segs, int *iseg, int nseg, int *itime, double T, double *xa)
{
    // Computes the next position, *xa, of a uniform sampled position profile.
    // The profile is composed of an array of segments (type: seg)
    // Each segment consists of:
    //      xfa:    final position
    //      v:      maximum velocity
    //      a:      maximum acceleration
    //      d:      dwell time at the final position
    //  Called from a loop, the profile proceeds from the current position,
    //  through each segment in turn, and then repeats.
    // Inputs:
    //  seg *segs:  - segments array
    //  int *iseg:  - variable hold segment index
    //  int nseg:   - number of segments in the profile
    //  int *itime  - time index within a segment (= -1 at segment beginning)
    //  double T:   - time increment
    // Outputs:
    //  double *xa: - next position in profile
    // Returns:     n - number of samples in the profile, 0 otherwise
	//
	//  Call with *itime = -1, *iseg = -1, outside the loop to initialize.

    double t, t1=0, t2=1, tf=1, tramp, x1=1, xramp, xfr=1, xr, d;
    static double x0, dir;
    static int ntot;
    double vmax=1, amax=1;
    int n;

    if (*itime==-1) {
        (*iseg)++;
        if(*iseg==nseg) {
        	*iseg=0;
        	ntot = 0;
        }
        *itime=0;
        x0=*xa;
    }
    vmax=segs[*iseg].v;
    amax=segs[*iseg].a;
    d=segs[*iseg].d;
    xfr=segs[*iseg].xfa-x0;
    dir=1.0;
    if(xfr<0){
        dir=-1.;
        xfr=-xfr;
    }
    t1 = vmax/amax;
    x1 = 1./2.*amax*t1*t1;
    if (x1<xfr/2) {
        xramp = xfr-2.*x1;
        tramp = xramp/vmax;
        t2 = t1+tramp;
        tf = t2+t1;
    } else {
        x1 = xfr/2;
        t1 = sqrt(2*x1/amax);
        t2 = t1;
        tf = 2.*t1;
    }
    n = trunc((tf+d)/T)+1;

    t = *itime*T;
    if(t<t1) {
        xr = 1./2.*amax*t*t;
    } else if (t>=t1 && t<t2) {
        xr = x1+vmax*(t-t1);
    } else if (t>=t2 && t<tf) {
        xr = xfr-1./2.*amax*(tf-t)*(tf-t);
    } else {
        xr = xfr;
    }
    *xa=x0+dir*xr;
    (*itime)++;
    if(*itime==n+1) {
    	ntot = ntot + *itime - 1;
        *itime=-1;
        if(*iseg==nseg-1) {
        	return ntot;
        }
    }
    return 0;
}
