#include"Graph.h"
#define _USE_MATH_DEFINES
#include<stdlib.h>
#include<math.h>
#include<string>
#include<glut.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<complex>

Graph::Graph()
{
	glLineStipple(1 , 0x0F0F);
	glLineWidth(3);

	offset = 0.1;
	X = Range<double>(100, 1000, 1);

	addLineData("");
}

Graph::~Graph(){

}


void Graph::GraphVertex(Complex a){

}

void Graph::GraphVertex(double x,double y){
	double r = 1.0 - offset;
	double wid = 2.0*r/(X.getMax() - X.getMin());
	double hei = 2.0*r/(Y.getMax() - Y.getMin());
	glVertex2d( wid * (x-X.getMin()) - r , hei * (y - Y.getMin() ) - r);
}

void Graph::addLineData(string str){
	string folder = "../../../../ë≤å§/SelfModelâêÕâÉvÉçÉOÉâÉÄ/d=50M=10n=(1.56,1)/";
	for(int i=0; i <= 45 ;i+= 15){
		string name = to_s(i) + "[deg].txt";
		DataSet<double> D = DataSet<double>(folder+name, X);
		Dash.push_back(D);
		Y.setMin( D.Y.getMin() );
		Y.setMax( D.Y.getMax() );
	}
}

void Graph::addSolidData(string str){
	string f = "ana/";
	string st[4] = {"400", "420", "450", "480"};

	string name = "";
	for(int i=0; i<4; i++){
		DataSet<double> D = DataSet<double>(name, X);
		Solid.push_back(D);
		Y.setMin( D.Y.getMin() );
		Y.setMax( D.Y.getMax() );
	}
}

void Graph::draw_data(DataSet<double> dat, bool f){
	glDisable(GL_LINE_STIPPLE);
	MyColor c = dat.color;

	glColor3d(c.r, c.g, c.b);

	if(f){
		glEnable(GL_LINE_STIPPLE);
	}
	glBegin(GL_LINE_STRIP);

	vector<complex<double>>::iterator it = dat.Data.begin();
	while(it != dat.Data.end()){
		GraphVertex((*it).real(), (*it).imag());
		it++;
	}

	glEnd();
}

void Graph::draw_axis(){
	glLineWidth(2.0);
	glColor3d(0.1,0.1,0.1);
	glDisable(GL_LINE_STIPPLE);

	glBegin(GL_LINES);
		GraphVertex( X.getMin(), Y.getMin() );
		GraphVertex( X.getMin(), Y.getMax() );
	glEnd();
	glBegin(GL_LINES);
		GraphVertex( X.getMax(), Y.getMin() );
		GraphVertex( X.getMax(), Y.getMax() );
	glEnd();
	
	glBegin(GL_LINES);
		GraphVertex(X.getMin(), Y.getMax());
		GraphVertex(X.getMax(), Y.getMax());
	glEnd();

	glBegin(GL_LINES);
		GraphVertex(X.getMin(), Y.getMin());
		GraphVertex(X.getMax(), Y.getMin());
	glEnd();
	glLineWidth(1.0);
}

void Graph::draw(){
//	glClear(GL_COLOR_BUFFER_BIT);
	glColor3d(10,10,10);		//îwåiÇîíÇ≈ï`âÊ
	glRectd(-1,-1, 1,1);

	draw_axis();
	glLineWidth(1.0);
	vector<DataSet<double>>::iterator it;
	
	it = Solid.begin();

	while(it != Solid.end() ){
		draw_data((*it), true);
		it++;
	}

	it = Dash.begin();
	while(it != Dash.end() ){
		draw_data((*it), false);
		it++;

	}

}