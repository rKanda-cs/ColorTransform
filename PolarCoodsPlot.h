#ifndef _POLAR_COODS_PLOT_H
#define _POLAR_COODS_PLOT_H

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include <string>
#include <glut.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "RGBTransform.h"
#include "Drawer.h"
using namespace std;

typedef vector<double> DataArray;

class PolarCoodsPlot: public Drawer{
	int dataNum;
	DataArray data;
	double range;
	double v_max, v_min;
public:
	PolarCoodsPlot(string name){

		ifstream fp(name);
		v_max = -100, v_min = 100;
		double tmp;
		dataNum = 0;
		while(fp){
			fp >> tmp;
			data.push_back(tmp);
			v_min = min(v_min, tmp);
			dataNum++;
		}

		for(int i=0; i<dataNum; i++){
			data[i] -= v_min;
			v_max = max(v_max, data[i]);
		}
		range = 0.9/v_max;

		cout << dataNum << endl;
		cout << v_max << endl;
		cout << v_min << endl;
	};

	void draw(){
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3d(255,255,255);		//”wŒi‚ð”’‚Å•`‰æ
		glRectd(-1,-1, 1,1);

		glLineWidth(1);
		glColor3d(0,0,0);

		//‰~‚Ì•`‰æ
		glBegin(GL_LINE_LOOP);
		for(int i=0; i < 360; i++){
			double rad = i/180.0*M_PI;
			glVertex2d(0.9*cos(rad) , 0.9*sin(rad));
		}
		glEnd();

		//ŠDF‚Å–Ú·‚è‚ð•`‰æ	DB•\‹L
		glColor3d(0.8, 0.8, 0.8);
		for(int i=0; i < 12; i++){
			double rad = i/6.0*M_PI;
			glBegin(GL_LINE_LOOP);
			glVertex2d(0,0);
			glVertex2d(0.9*cos(rad) , 0.9*sin(rad));
			glEnd();
		}
		/*
		//Ô‚Åƒf[ƒ^‚ð•`‰æ
		glColor3d(255,0,0);
		for(int i=0; i < 5; i++){
			double v = i*v_max/5;
			double s = v*range;

			drawBitmapString(GLUT_BITMAP_HELVETICA_18, to_s(s+v_min), 0, s);
		}
		*/
		glColor3d(255,0,0);
		glBegin(GL_LINE_LOOP);
		for(int i=0; i < dataNum; i++){
			double rad = i/180.0*M_PI;
			double x = range*(data[i]*cos(rad));
			double y = range*(data[i]*sin(rad));
			glVertex2d(x , y);
		}
		glEnd();
	};

};

#endif //_POLAR_COODS_PLOT_H