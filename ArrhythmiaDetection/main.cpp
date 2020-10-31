//
//  main.cpp
//  ArrhythmiaDetection
//
//  Created by Abdulrahman on 4/4/15.
//  Copyright (c) 2015 Abdulrahman. All rights reserved.
//

#include  <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "qrsdet.h"
#include <stdlib.h>
#include <time.h>

using namespace std;

int QRSFilter(int datum,int init);
int Peak( int datum, int init );

void checkTachycardiaBradycardia();
void RDetection();
void timeSignalWriter();

//Added
#include <math.h>
#include <cstring>

#include "ArrhythmiaClassification.h"

#define PRE_BLANK	MS200

int deriv1(int x0, int init);


// Local Prototypes.

int Peak(int datum, int init);
int median(int* array, int datnum);
int thresh(int qmedian, int nmedian);
int BLSCheck(int* dBuf, int dbPtr, int* maxder);
void RDetection();

int lpfilt(int datum, int init);
int hpfilt(int datum, int init);
int deriv1(int x0, int init);
int deriv2(int x0, int init);
int mvwint(int datum, int init);


double TH = 0.475;

int DDBuffer[DER_DELAY], DDPtr;	/* Buffer holding derivative data. */
int Dly = 0;

const int MEMMOVELEN = 7 * sizeof(int);

int QRSDet(int datum, int init)
{
	//Added
	cout << "QRSDet function runs.....\n";

	static int det_thresh, qpkcnt = 0;
	static int qrsbuf[8], noise[8], rrbuf[8];
	static int rsetBuff[8], rsetCount = 0;
	static int nmedian, qmedian, rrmedian;
	static int count, sbpeak = 0, sbloc, sbcount = MS1500;
	static int maxder, lastmax;
	static int initBlank, initMax;
	static int preBlankCnt, tempPeak;

	int fdatum, QrsDelay = 0;
	int i, newPeak, aPeak;

	/*	Initialize all buffers to 0 on the first call.	*/

	if (init)
	{
		for (i = 0; i < 8; ++i)
		{
			noise[i] = 0;	/* Initialize noise buffer */
			rrbuf[i] = MS1000;/* and R-to-R interval buffer. */
		}
			
		qpkcnt = maxder = lastmax = count = sbpeak = 0;
		initBlank = initMax = preBlankCnt = DDPtr = 0;
		sbcount = MS1500;
		QRSFilter(0, 1);	/* initialize filters. */
		Peak(0, 1);
	}

	fdatum = QRSFilter(datum, 0);	/* Filter data. */


	/* Wait until normal detector is ready before calling early detections. */

	aPeak = Peak(fdatum, 0);

	// Hold any peak that is detected for 200 ms
	// in case a bigger one comes along.  There
	// can only be one QRS complex in any 200 ms window.

	newPeak = 0;
	if (aPeak && !preBlankCnt)			// If there has been no peak for 200 ms
	{										// save this one and start counting.
		tempPeak = aPeak;
		preBlankCnt = PRE_BLANK;			// MS200
	}

	else if (!aPeak && preBlankCnt)	// If we have held onto a peak for
	{										// 200 ms pass it on for evaluation.
		if (--preBlankCnt == 0)
			newPeak = tempPeak;
	}

	else if (aPeak)							// If we were holding a peak, but
	{										// this ones bigger, save it and
		if (aPeak > tempPeak)				// start counting to 200 ms again.
		{
			tempPeak = aPeak;
			preBlankCnt = PRE_BLANK; // MS200
		}
		else if (--preBlankCnt == 0)
			newPeak = tempPeak;
	}

	/*	newPeak = 0 ;
		if((aPeak != 0) && (preBlankCnt == 0))
			newPeak = aPeak ;
		else if(preBlankCnt != 0) --preBlankCnt ; */



		/* Save derivative of raw signal for T-wave and baseline
		   shift discrimination. */

	DDBuffer[DDPtr] = deriv1(datum, 0);
	if (++DDPtr == DER_DELAY)
		DDPtr = 0;

	/* Initialize the qrs peak buffer with the first eight 	*/
	/* local maximum peaks detected.						*/

	if (qpkcnt < 8)
	{
		++count;
		if (newPeak > 0) count = WINDOW_WIDTH;
		if (++initBlank == MS1000)
		{
			initBlank = 0;
			qrsbuf[qpkcnt] = initMax;
			initMax = 0;
			++qpkcnt;
			if (qpkcnt == 8)
			{
				qmedian = median(qrsbuf, 8);
				nmedian = 0;
				rrmedian = MS1000;
				sbcount = MS1500 + MS150;
				det_thresh = thresh(qmedian, nmedian);
			}
		}
		if (newPeak > initMax)
			initMax = newPeak;
	}

	else	/* Else test for a qrs. */
	{
		++count;
		if (newPeak > 0)
		{


			/* Check for maximum derivative and matching minima and maxima
			   for T-wave and baseline shift rejection.  Only consider this
			   peak if it doesn't seem to be a base line shift. */

			if (!BLSCheck(DDBuffer, DDPtr, &maxder))
			{


				// Classify the beat as a QRS complex
				// if the peak is larger than the detection threshold.

				if (newPeak > det_thresh)
				{
					memmove(&qrsbuf[1], qrsbuf, MEMMOVELEN);
					qrsbuf[0] = newPeak;
					qmedian = median(qrsbuf, 8);
					det_thresh = thresh(qmedian, nmedian);
					memmove(&rrbuf[1], rrbuf, MEMMOVELEN);
					rrbuf[0] = count - WINDOW_WIDTH;
					rrmedian = median(rrbuf, 8);
					sbcount = rrmedian + (rrmedian >> 1) + WINDOW_WIDTH;
					count = WINDOW_WIDTH;

					sbpeak = 0;

					lastmax = maxder;
					maxder = 0;
					QrsDelay = WINDOW_WIDTH + FILTER_DELAY;
					initBlank = initMax = rsetCount = 0;

					//		preBlankCnt = PRE_BLANK ;
				}

				// If a peak isn't a QRS update noise buffer and estimate.
				// Store the peak for possible search back.


				else
				{
					memmove(&noise[1], noise, MEMMOVELEN);
					noise[0] = newPeak;
					nmedian = median(noise, 8);
					det_thresh = thresh(qmedian, nmedian);

					// Don't include early peaks (which might be T-waves)
					// in the search back process.  A T-wave can mask
					// a small following QRS.

					if ((newPeak > sbpeak) && ((count - WINDOW_WIDTH) >= MS360))
					{
						sbpeak = newPeak;
						sbloc = count - WINDOW_WIDTH;
					}
				}
			}
		}

		/* Test for search back condition.  If a QRS is found in  */
		/* search back update the QRS buffer and det_thresh.      */

		if ((count > sbcount) && (sbpeak > (det_thresh >> 1)))
		{
			memmove(&qrsbuf[1], qrsbuf, MEMMOVELEN);
			qrsbuf[0] = sbpeak;
			qmedian = median(qrsbuf, 8);
			det_thresh = thresh(qmedian, nmedian);
			memmove(&rrbuf[1], rrbuf, MEMMOVELEN);
			rrbuf[0] = sbloc;
			rrmedian = median(rrbuf, 8);
			sbcount = rrmedian + (rrmedian >> 1) + WINDOW_WIDTH;
			QrsDelay = count = count - sbloc;
			QrsDelay += FILTER_DELAY;
			sbpeak = 0;
			lastmax = maxder;
			maxder = 0;
			initBlank = initMax = rsetCount = 0;
		}
	}

	// In the background estimate threshold to replace adaptive threshold
	// if eight seconds elapses without a QRS detection.

	if (qpkcnt == 8)
	{
		if (++initBlank == MS1000)
		{
			initBlank = 0;
			rsetBuff[rsetCount] = initMax;
			initMax = 0;
			++rsetCount;

			// Reset threshold if it has been 8 seconds without
			// a detection.

			if (rsetCount == 8)
			{
				for (i = 0; i < 8; ++i)
				{
					qrsbuf[i] = rsetBuff[i];
					noise[i] = 0;
				}
				qmedian = median(rsetBuff, 8);
				nmedian = 0;
				rrmedian = MS1000;
				sbcount = MS1500 + MS150;
				det_thresh = thresh(qmedian, nmedian);
				initBlank = initMax = rsetCount = 0;
				sbpeak = 0;
			}
		}
		if (newPeak > initMax)
			initMax = newPeak;
	}

	return(QrsDelay);
}

/**************************************************************
* peak() takes a datum as input and returns a peak height
* when the signal returns to half its peak height, or
**************************************************************/

int Peak(int datum, int init)
{
	//Added
	cout << "Peak function runs.....\n";

	static int max = 0, timeSinceMax = 0, lastDatum;
	int pk = 0;

	if (init)
		max = timeSinceMax = 0;

	if (timeSinceMax > 0)
		++timeSinceMax;

	if ((datum > lastDatum) && (datum > max))
	{
		max = datum;
		if (max > 2)
			timeSinceMax = 1;
	}

	else if (datum < (max >> 1))
	{
		pk = max;
		max = 0;
		timeSinceMax = 0;
		Dly = 0;
	}

	else if (timeSinceMax > MS95)
	{
		pk = max;
		max = 0;
		timeSinceMax = 0;
		Dly = 3;
	}
	lastDatum = datum;
	return(pk);
}

/********************************************************************
median returns the median of an array of integers.  It uses a slow
sort algorithm, but these arrays are small, so it hardly matters.
********************************************************************/

int median(int* array, int datnum)
{
	//Added
	cout << "median function runs.....\n";

	int i, j, k, temp, sort[20];
	for (i = 0; i < datnum; ++i)
		sort[i] = array[i];
	for (i = 0; i < datnum; ++i)
	{
		temp = sort[i];
		for (j = 0; (temp < sort[j]) && (j < i); ++j);
		for (k = i - 1; k >= j; --k)
			sort[k + 1] = sort[k];
		sort[j] = temp;
	}
	return(sort[datnum >> 1]);
}
/*
int median(int *array, int datnum)
	{
	long sum ;
	int i ;

	for(i = 0, sum = 0; i < datnum; ++i)
		sum += array[i] ;
	sum /= datnum ;
	return(sum) ;
	} */

	/****************************************************************************
	 thresh() calculates the detection threshold from the qrs median and noise
	 median estimates.
	****************************************************************************/

int thresh(int qmedian, int nmedian)
{
	//Added
	cout << "thresh function runs.....\n";

	int thrsh, dmed;
	double temp;
	dmed = qmedian - nmedian;
	/*	thrsh = nmedian + (dmed>>2) + (dmed>>3) + (dmed>>4); */
	temp = dmed;
	temp *= TH;
	dmed = temp;
	thrsh = nmedian + dmed; /* dmed * THRESHOLD */
	return(thrsh);
}

/***********************************************************************
	BLSCheck() reviews data to see if a baseline shift has occurred.
	This is done by looking for both positive and negative slopes of
	roughly the same magnitude in a 220 ms window.
***********************************************************************/

int BLSCheck(int* dBuf, int dbPtr, int* maxder)
{
	//Added
	cout << "BLSCheck function runs.....\n";

	int max, min, maxt, mint, t, x;
	max = min = 0;

	return(0);

	for (t = 0; t < MS220; ++t)
	{
		x = dBuf[dbPtr];
		if (x > max)
		{
			maxt = t;
			max = x;
		}
		else if (x < min)
		{
			mint = t;
			min = x;
		}
		if (++dbPtr == DER_DELAY)
			dbPtr = 0;
	}

	*maxder = max;
	min = -min;

	/* Possible beat if a maximum and minimum pair are found
		where the interval between them is less than 150 ms. */

	if ((max > (min >> 3)) && (min > (max >> 3)) &&
		(fabs(maxt - mint) < MS150))
		return(0);

	else
		return(1);
}
/********************************************************************
 RDetection() a simple function used for Writing R value of QRS
 complex.
 ********************************************************************/
void RDetection() {

	//Added
	cout << "RDetection function runs.....\n";

	ofstream RFile("Rfile.txt");
	ifstream QrsDelay("QRSDelay.txt");
	bool qrsComplex = false;
	bool rFound = false;

	int id, signal, qValue;
	double time;

	int previos = 0, currentSignal = 0, init = 0;
	int previosID = 0, previosSignal = 0, previosQvalue = 0;

	while (QrsDelay >> id >> time >> signal >> qValue) {
		if (init) {
			if (previosQvalue)qrsComplex = true;
			if (qrsComplex) {
				if (currentSignal > previos && currentSignal > signal) {
					qrsComplex = false;
					RFile << to_string(previosID)
						<< "\t" << to_string(previosSignal)
						<< "\t" << to_string(previosQvalue)
						<< "\t" << to_string(currentSignal) << "\n";
					currentSignal = 0;
					previos = 0;
					previosQvalue = qValue;
					previosSignal = signal;
					previosID = id;

				}
				else {
					previos = currentSignal;
					currentSignal = signal;

					RFile << to_string(previosID)
						<< "\t" << to_string(previosSignal)
						<< "\t" << to_string(previosQvalue)
						<< "\t0\n";
					previosQvalue = qValue;
					previosSignal = signal;
					previosID = id;
				}

			}
			else {
				RFile << to_string(previosID)
					<< "\t" << to_string(previosSignal)
					<< "\t" << to_string(previosQvalue)
					<< "\t0\n";
				previosQvalue = qValue;
				previosSignal = signal;
				previosID = id;
			}
		}
		else {/*if initilaizer*/
			previosID = id;
			previosQvalue = qValue;
			previosSignal = signal;
			init = 1;
		}
	}
}

/******************************************************************************
* Syntax:
*	int QRSFilter(int datum, int init) ;
* Description:
*	QRSFilter() takes samples of an ECG signal as input and returns a sample of
*	a signal that is an estimate of the local energy in the QRS bandwidth.  In
*	other words, the signal has a lump in it whenever a QRS complex, or QRS
*	complex like artifact occurs.  The filters were originally designed for data
*  sampled at 200 samples per second, but they work nearly as well at sample
*	frequencies from 150 to 250 samples per second.
*
*	The filter buffers and static variables are reset if a value other than
*	0 is passed to QRSFilter through init.
*******************************************************************************/

int QRSFilter(int datum, int init)
{
	//Added
	cout << "QRSFilter function runs.....\n";

	int fdatum;

	if (init)
	{
		hpfilt(0, 1);		// Initialize filters.
		lpfilt(0, 1);
		mvwint(0, 1);
		deriv1(0, 1);
		deriv2(0, 1);
	}

	fdatum = lpfilt(datum, 0);		// Low pass filter data.
	fdatum = hpfilt(fdatum, 0);	// High pass filter data.
	fdatum = deriv2(fdatum, 0);	// Take the derivative.
	fdatum = fabs(fdatum);				// Take the absolute value.
	//fdatum = mvwint( fdatum, 0 ) ;	// Average over an 80 ms window .
	return(fdatum);
}


/*************************************************************************
*  lpfilt() implements the digital filter represented by the difference
*  equation:
*
* 	y[n] = 2*y[n-1] - y[n-2] + x[n] - 2*x[t-24 ms] + x[t-48 ms]
*
*	Note that the filter delay is (LPBUFFER_LGTH/2)-1
*
**************************************************************************/

int lpfilt(int datum, int init)
{
	//Added
	cout << "lpfilt function runs.....\n";

	static long y1 = 0, y2 = 0;
	static int data[LPBUFFER_LGTH], ptr = 0;
	long y0;
	int output, halfPtr;
	if (init)
	{
		for (ptr = 0; ptr < LPBUFFER_LGTH; ++ptr)
			data[ptr] = 0;
		y1 = y2 = 0;
		ptr = 0;
	}
	halfPtr = ptr - (LPBUFFER_LGTH / 2);	// Use halfPtr to index
	if (halfPtr < 0)							// to x[n-6].
		halfPtr += LPBUFFER_LGTH;
	y0 = (y1 << 1) - y2 + datum - (data[halfPtr] << 1) + data[ptr];
	y2 = y1;
	y1 = y0;
	output = y0 / ((LPBUFFER_LGTH * LPBUFFER_LGTH) / 4);
	data[ptr] = datum;			// Stick most recent sample into
	if (++ptr == LPBUFFER_LGTH)	// the circular buffer and update
		ptr = 0;					// the buffer pointer.
	return(output);
}


/******************************************************************************
*  hpfilt() implements the high pass filter represented by the following
*  difference equation:
*
*	y[n] = y[n-1] + x[n] - x[n-128 ms]
*	z[n] = x[n-64 ms] - y[n] ;
*
*  Filter delay is (HPBUFFER_LGTH-1)/2
******************************************************************************/

int hpfilt(int datum, int init)
{
	//Added
	cout << "hpfilt function runs.....\n";

	static long y = 0;
	static int data[HPBUFFER_LGTH], ptr = 0;
	int z, halfPtr;

	if (init)
	{
		for (ptr = 0; ptr < HPBUFFER_LGTH; ++ptr)
			data[ptr] = 0;
		ptr = 0;
		y = 0;
	}

	y += datum - data[ptr];
	halfPtr = ptr - (HPBUFFER_LGTH / 2);
	if (halfPtr < 0)
		halfPtr += HPBUFFER_LGTH;
	z = data[halfPtr] - (y / HPBUFFER_LGTH);

	data[ptr] = datum;
	if (++ptr == HPBUFFER_LGTH)
		ptr = 0;

	return(z);
}

/*****************************************************************************
*  deriv1 and deriv2 implement derivative approximations represented by
*  the difference equation:
*
*	y[n] = x[n] - x[n - 10ms]
*
*  Filter delay is DERIV_LENGTH/2
*****************************************************************************/

int deriv1(int x, int init)
{
	//Added
	cout << "deriv1 function runs.....\n";

	static int derBuff[DERIV_LENGTH], derI = 0;
	int y;

	if (init != 0)
	{
		for (derI = 0; derI < DERIV_LENGTH; ++derI)
			derBuff[derI] = 0;
		derI = 0;
		return(0);
	}

	y = x - derBuff[derI];
	derBuff[derI] = x;
	if (++derI == DERIV_LENGTH)
		derI = 0;
	return(y);
}

int deriv2(int x, int init)
{
	//Added
	cout << "deriv2 function runs.....\n";

	static int derBuff[DERIV_LENGTH], derI = 0;
	int y;

	if (init != 0)
	{
		for (derI = 0; derI < DERIV_LENGTH; ++derI)
			derBuff[derI] = 0;
		derI = 0;
		return(0);
	}

	y = x - derBuff[derI];
	derBuff[derI] = x;
	if (++derI == DERIV_LENGTH)
		derI = 0;
	return(y);
}




/*****************************************************************************
* mvwint() implements a moving window integrator.  Actually, mvwint() averages
* the signal values over the last WINDOW_WIDTH samples.
*****************************************************************************/

int mvwint(int datum, int init)
{
	//Added
	cout << "mvwint function runs.....\n";

	static long sum = 0;
	static int data[WINDOW_WIDTH], ptr = 0;
	int output;
	if (init)
	{
		for (ptr = 0; ptr < WINDOW_WIDTH; ++ptr)
			data[ptr] = 0;
		sum = 0;
		ptr = 0;
	}
	sum += datum;
	sum -= data[ptr];
	data[ptr] = datum;
	if (++ptr == WINDOW_WIDTH)
		ptr = 0;
	if ((sum / WINDOW_WIDTH) > 32000)
		output = 32000;
	else
		output = sum / WINDOW_WIDTH;
	return(output);
}

//
//  fileHandling.cpp
//  ArrhythmiaDetection
//
//  Created by Abdulrahman on 4/5/15.
//  Copyright (c) 2015 Abdulrahman. All rights reserved.
//

void timeSignalWriter() {

	//Added
	cout << "timeSignalWriter function runs.....\n";

	ifstream signal("standData.txt");
	ifstream csvBased("csvData.txt");
	ofstream timedSignal("timedSignal.txt");



	double time, lowSignal, highSignal;
	int id, high, low;
	while ((signal >> id >> high >> low) && (csvBased >> time >> lowSignal >> highSignal)) {


		string idS = to_string(id);
		string timeS = to_string(time);
		string highS = to_string(high);
		string highS2 = to_string(highSignal);
		string lowS = to_string(low);

		timedSignal << idS << "\t" << timeS << "\t" << highS << "\n";


	}
	timedSignal.close();
}

//
//  ArrhythmiaClassification.cpp
//  ArrhythmiaDetection
//
//  Created by Abdulrahman on 4/4/15.
//  Copyright (c) 2015 Abdulrahman. All rights reserved.
//


int TA = 0, NR = 0, BR = 0;
float startTime = 0;
void checkTachycardiaBradycardia(int peak, float time) {

	//Added
	cout << "checkTachycardiaBradycardia runs......." << endl;
	

	if (peak) {
		double heartRate = 60 / (time - startTime);
		printf("\t\t\t\t\t\n\nHeart rate : %f\n\n", heartRate);
		if (heartRate > 100) {
			TA++;
		}
		else if (heartRate < 60) {
			BR++;
		}
		else
		{	
			NR++;
		}
		startTime = time;
	}
}




int main() {

    //Added
    cout << "Main function runs.....\n";

    ofstream signal("QRSDelay.txt");
    ofstream timerReport("timeReport.txt");

    timeSignalWriter();


    ifstream timedSignal("timedSignal.txt");

    bool start = false;
    double a, c;
    double b;
    int init = 1, i = 0;
    //double startTime = 0.0 , endTime = 0.0;
    string s;
    int signals[1000];
    int QRSPeak[1000];
    int ids[1000];
    double times[1000];


    const clock_t timer = clock();
    while (timedSignal >> a >> b >> c) {

        if (init) {
            QRSDet(0, 1);
            init = 0;
        }

        int qrsDelay = QRSDet(c, init);

		//printf("QrsDelay : %lf, %lf, %lf, %lf\n", qrsDelay, a, b, c);

        if (qrsDelay) {
            start = true;
        }

        if (start) {
            if (qrsDelay) {

                QRSPeak[i - qrsDelay] = signals[i - qrsDelay];
                for (int j = 0; j < i; j++) {
                    s = to_string(ids[j]);
                    signal << s << "\t";
                    s = to_string(times[j]);
                    signal << s << "\t";
                    s = to_string(signals[j]);
                    signal << s << "\t";
                    s = to_string(QRSPeak[j]);
                    signal << s << "\n";

					checkTachycardiaBradycardia(QRSPeak[j], times[j]);
                }
                i = 0;

                //Added
                memset(QRSPeak, 0, sizeof(QRSPeak));
                memset(signals, 0, sizeof(signals));
                memset(times, 0, sizeof(times));
                memset(ids, 0, sizeof(ids));

            }
            else if (i == 1000) {
                for (int j = 0; j < 500; j++) {
                    s = to_string(ids[j]);
                    ids[j] = ids[j + 500];
                    signal << s << "\t";
                    s = to_string(times[j]);
					
					float t = times[j];

                    times[j] = times[j + 500];
                    signal << s << "\t";
                    s = to_string(signals[j]);
                    signals[j] = signals[j + 500];
                    signal << s << "\t";
                    s = to_string(QRSPeak[j]);

					checkTachycardiaBradycardia(QRSPeak[j], t);
                    
					QRSPeak[j] = QRSPeak[j + 500];
                    signal << s << "\n";
                }
                i = 500;
            }
            else {
                ids[i] = a;
                times[i] = b;
                signals[i] = c;
                QRSPeak[i] = qrsDelay;
                i++;

            }

        }
    }
    float detectionTime = float((clock() - timer) / CLOCKS_PER_SEC);
    /*checkTachycardiaBradycardia();*/
    float arrthmieaDTime = float((clock() - timer) / CLOCKS_PER_SEC) - detectionTime;
    timerReport << "Timer Report For Record [103] :-\n";
    timerReport << "QRS Detection Time =\t" << to_string(detectionTime) << "s\n";
    timerReport << "ARR Detection Time =\t" << to_string(arrthmieaDTime) << "s\n";
    timerReport << "Total Time =\t" << to_string(timer) << "s\n";

    signal.close();
    timerReport.close();
    RDetection();

    //Added
	printf("TA : %d\tBR : %d\tNR : %d \n", TA, BR, NR);

    cout << "\n\n\t\t\t\tMain function ends.....\n\n\n\n\n\n";

    return 0;
}
