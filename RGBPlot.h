#ifndef _RGBPLOT_H
#define _RGBPLOT_H
#include"Drawer.h"
#include"RGBTransform.h"
#include<string>
#include<vector>
#include<map>
#include<glut.h>
#include<complex>

#define PI 3.141592

const int out_min_deg = 0;	//最小角度
const int out_max_deg = 180;	//最大角度
const int out_d_deg   = 1;		//角度間隔
const int out_num	  = (out_max_deg-out_min_deg) / out_d_deg;
const int in_min_deg  = 0;
const int in_max_deg  = 180;
const int in_d_deg    = 5;
const int in_num	  = (in_max_deg-in_min_deg) / in_d_deg;
const int cellNum	  = (in_num+1)*(out_num+1);
const double offset = 0.2;

bool between(int mi, int ma, int a){
	return  (a>=mi && a<=ma);
}

class RGBPlot:public Drawer{
RGBTransform *rgb;
MyColor *ColorMap;
double max_sum;
string Dir;
public:
	RGBPlot() {
		rgb = new RGBTransform();
		ColorMap = new MyColor[cellNum];

		for (int i = 0; i < cellNum; i++)
			ColorMap[i] = MyColor(0.0, 0.0, 0.0);

		max_sum = 0;

		setDirTM();
		ifstream fp(getDir() + "WaveAngleStrength.txt");				//ここから
		while (fp) {
			double tmp;
			complex<double> tmp1;
			fp >> tmp1 >> tmp;
			max_sum = max(max_sum, tmp);	//最大の総和で正規化
		}																//ここまで　matlab理論値の時いらない

		/*
		setDirTE();
		fp.open(getDir() + "WaveAngleStrength.txt" );
		while(fp){
			double tmp;
			complex<double> tmp1;
			fp >> tmp1 >> tmp;
			max_sum = max(max_sum, tmp);	//最大の総和で正規化
		}

		for(int j=-90; j <= 0; j+=10)
			for(int i=380; i<=700; i+=5)
				getPlot(i,j);

		setDirTM();
		*/
		//max_sum = 1;
		for(int j=0; j <= 180; j+=5){
			for(int i=380; i<=700; i+=5){
				//if((i>=600 && i<=700)) continue;
				getPlot(i,j);
			}
		}

//		getPlot_hair();
/*
		for (int i = 380; i <= 700; i++) {
			getPlot_theo(i);
		}
*/
//		saveColorData1();
//		saveColorData2();
//		saveColorData3();
		saveColorData4();

	}
	
	string getDir(){
		return Dir;
	}

	void setDirTE(){
		Dir = "../../../DataSet/TE/Morpho(1,1.56)M=8/60nm(nonShelf)(10nm,200cell)/NTFF/";	//20[nm]：緑, 50[nm]:緑
		//Dir = "../../../DataSet/TE/Morpho/90nm(nonShelf)(10nm,200cell)/NTFF/";
	}

	void setDirTM(){
		//Dir = "../../TheoreticalValue/data/";
		//Dir = "../../FDTD_HairSimulation/DataSet/HairModel/incidenceLayer_withSig/(50nm,640cell)/TM/Ns/NTFF/";
		//Dir = "../../FDTD_HairSimulation/DataSet/HairModel/incidenceLayer/(50nm,960cell)/TM/Ns/NTFF/";
		//Dir = "../../FDTD_HairSimulation/DataSet/HairModel/normalplane/e=0/(100nm,1280cell)/TM/Ns/NTFF/";
		//Dir = "../../FDTD_HairSimulation/DataSet/SlabModel/(10nm,586cell)/TM/Ns/7枚slab②_2018.07.01/NTFF/";
		Dir = "../../FDTD_HairSimulation/DataSet/SlabModel/(20nm,375cell)/TM/Ns/NTFF/";
		//Dir = "../../FDTD_HairSimulation/DataSet/SlabModel/scratch/(20nm,500cell)/TM/Ns/2019.1.9_3層12,6deg/NTFF/";
		//Dir = "../../FDTD_HairSimulation/DataSet/MorphoModel/(10nm,400cell)/TM/Ns/NTFF/";
		//Dir = "../../FDTD_HairSimulation/DataSet/SlabModel/(10nm,100cell)/TM/Ns/NTFF/";
		//Dir = "../../FDTD-Cpp-master/DataSet/Morpho(1,1.56)M=8/120nm(nonShelf)(10nm,200cell)/TM/St/NTFF/";
		//Dir = "../../../DataSet/ShelfModel/d=235M=9/";
		//Dir = "../../../DataSet/TM/Morpho(1.56,1)M=8/110nm(nonShelf)(10nm,200cell)/NTFF/";	//20[nm]：緑, 50[nm]:緑
		//Dir = "../../../DataSet/TM/Morpho/90nm(nonShelf)(10nm,200cell)/NTFF/";
	}
	void saveColorData(){
		ofstream ofp("Colordata.txt");
		for(int i=0; i< cellNum ;i++)
			ofp << ColorMap[i].r << "  "<< ColorMap[i].g << "  " << ColorMap[i].b << endl;
	}

	void saveColorData1() {						// 1次元カラーデータ
		ofstream ofp("Colordata1.txt");
		for (int i = 0; i < in_num; i++)
			ofp << Color(i, 0).r << "\n" << Color(i, 0).g << "\n" << Color(i, 0).b << endl;
		cout << "output Colordata1.txt" << endl;
	}

	void saveColorData2() {						// 2次元カラーデータ
		ofstream ofp("Colordata2.txt");
		int count = 0;
		for (int i = 0; i < in_num; i++) {
			for (int j = 0; j < out_num; j++) {
				ofp << Color(i, j).r << " " << Color(i, j).g << " " << Color(i, j).b << " " << i*in_d_deg << " " << j*out_d_deg << " " << count << endl;
				count++;
			}
		}
		cout << "output Colordata2.txt" << endl;
	}

	void saveColorData3() {						// 2次元カラーデータ	入射角反射角0-90度に整列
		ofstream ofp("Colordata3.txt");
		int count = 0;
		for (int i = in_num / 2; i <= in_num; i++) {
			for (int j = out_num / 2; j >= 0; j--) {
				ofp << Color(i, j).r << " " << Color(i, j).g << " " << Color(i, j).b << " " << (i*in_d_deg) - (in_num/2*in_d_deg) << " " << -((j*out_d_deg) - (out_num/2*out_d_deg)) << " " << count << endl;
				count++;
			}
		}
		cout << "output Colordata3.txt" << endl;
	}

	void saveColorData4() {						// 2次元カラーデータ	データがない空間を補間（出力画面と同じ）
		ofstream ofp("Colordata.txt");
		int count = 0;
		for (int i = 0; i < in_num; i++) {
			for (double k = 0.0; k < in_d_deg; k += 1) {
				for (int j = 0; j < out_num; j++) {
					double r = 1.0 / in_d_deg;
					double d = r*(k + 0.5);
					MyColor c1 = (1.0 - d) * Color(i, j) + d * Color(i + 1, j);
					MyColor c2 = (1.0 - d) * Color(i, j + 1) + d * Color(i + 1, j + 1);
					MyColor c = 0.5*(c1 + c2);

					ofp << c.r << " " << c.g << " " << c.b << " " << i*in_d_deg + k << " " << j*out_d_deg << " " << count << endl;
					count++;
				}
			}
		}
		cout << "output Colordata4.txt" << endl;
	}

	void getPlot(int lam, int in_deg){
		//string name   = to_s(in_deg) + "deg" + to_s(lam) + "nm.txt"; 
		string name = "(WL=" + to_s(lam) + "nm,AOI=" + to_s(in_deg) + "deg).txt"; 
		//string name = "(WL=" + to_s(lam) + "nm,AOI=" + to_s(135) + "deg).txt";

		ifstream fp(getDir()+name);	//ファイルを開く
		if(!fp)	cout <<  "file error"  << name << endl;

		if(in_deg < in_min_deg) return;
		int in = (in_deg - in_min_deg)/in_d_deg;
		int out_deg=0;	//散乱波の角度
		double ref;		//反射率

		while(fp){
			fp >> ref;
			ref = ref / max_sum;
			if( (out_deg < out_max_deg) && (out_deg >= out_min_deg)){
				double out = (out_deg - out_min_deg)/out_d_deg;
				Color(in,out) = Color(in,out) + rgb->CIE_RGB(rgb->getXYZ(lam, 3*ref));
			}
			out_deg++;
		}
	}

	void getPlot_theo(int lam) {
		string name = "TheoreticalVal_" + to_s(lam) + "nm.txt";

		ifstream fp(getDir() + name);
		if (!fp)	cout << "file error" << name << endl;

		int in = 0;
		double ref;
		double sum = 0;
		
		fp >> sum;
		while (fp) {		// in = out = 0～90 / 1deg
			fp >> ref;
			ref = ref / sum;
			for (int out = 0; out <= out_max_deg; out++) {
				Color(in, out) = Color(in, out) + rgb->CIE_RGB(rgb->getXYZ(lam, ref));
			}
			in++;
		}
	}

	void getPlot_hair() {
		MyColor *original = new MyColor[cellNum];
		for (int i = 0; i < cellNum; i++)
			original[i] = ColorMap[i];

		double in_angle = 90;		// 照射角度
		double r = out_num / 2;

		for (int j = 0; j <= out_num; j++) {
			double y = j - 90;
			double theta = asin(y / r);
			theta = theta * 180 / PI;

			int l_in_angle = in_angle - theta + 90;
			int l_out_angle = -theta + 90;

			int inc = 5;			// キューティクルの傾きを考慮
			l_in_angle -= inc;
			l_out_angle -= inc;

			if (l_in_angle >= 180) {	// 光が当たらない面
				for (int i = 0; i <= in_num; i++) {
					int a = (out_num + 1)*i + j;
					ColorMap[a] = MyColor(0, 0, 0);
				}
			}
			else if(l_in_angle < 180) {
				int x = (out_num + 1)*(l_in_angle / in_d_deg) + (l_out_angle / out_d_deg);
				for (int i = 0; i <= in_num; i++) {
					int a = (out_num + 1)*i + j;
					ColorMap[a] = original[x];
				}
			}
		}
	}

	
	void draw_axis(){
		double dh = 1.0-offset;
		double wid = 2.0*dh/in_num;
		double hei = 2.0*dh/out_num;
		glColor3d(0.0, 0.0, 0.0);

		// x軸ラベル

//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, " 0", -dh-0.15*offset , -dh-0.4*offset);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, " 180",  dh-0.15*offset , -dh-0.4*offset);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "theta(deg)", -0.1, -dh -0.1);

		drawBitmapString(GLUT_BITMAP_HELVETICA_12, " -90", -dh - 0.15*offset, -dh - 0.4*offset);
		drawBitmapString(GLUT_BITMAP_HELVETICA_12, " 0",  - 0.15*offset, -dh - 0.4*offset);
		drawBitmapString(GLUT_BITMAP_HELVETICA_12, " 90", dh - 0.15*offset, -dh - 0.4*offset);
		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "theta(deg)", -0.1, -dh - 0.15);

//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "incidence 135deg", -0.1, -dh - 0.1);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "light position = 90deg", -0.1, -dh - 0.1);


		//y軸ラベル

/*		for (int i = out_min_deg; i <= out_max_deg; i += 30) {
			double y = (i - out_min_deg) / out_d_deg*hei - dh;
			drawBitmapString(GLUT_BITMAP_HELVETICA_12, to_s(i - 270), -1.0 + 0.5*offset, y);
		}
*/
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "0", -dh - 0.5*offset,-dh-0.15*offset);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "45", -dh - 0.5*offset, 45 * hei - dh);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "135", -dh - 0.5*offset, 135 * hei - dh);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "180", -dh - 0.5*offset, dh-0.15*offset);
//		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "phi (deg)"  , -dh - 0.15, -0.01);

		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "-90", -dh - 0.5*offset, -dh - 0.15*offset);
		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "0", -dh - 0.5*offset, -0.15*offset);
		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "90", -dh - 0.5*offset, dh - 0.15*offset);
		drawBitmapString(GLUT_BITMAP_HELVETICA_12, "phi (deg)", -dh - 0.15, -0.1);

/*
		//外枠線(黒)
		glColor3d(0.0, 0.0, 0.0);
		glBegin(GL_LINES);
		glVertex2d(-dh, -dh);	glVertex2d( dh, -dh);
		glVertex2d(-dh, -dh);	glVertex2d(-dh,  dh);
		glVertex2d(-dh,  dh);	glVertex2d( dh,  dh);
		glVertex2d( dh, -dh);	glVertex2d( dh,  dh);
		glEnd();
*/
		// 白線目盛り
		glColor3d(1.0, 1.0, 1.0);
		glBegin(GL_LINES);
		for(int i=-1; i <= 1; i++){
			glVertex2d(-dh, 0.5*i*dh);
			glVertex2d( dh, 0.5*i*dh);
			glVertex2d(0.5*i*dh, -dh);
			glVertex2d(0.5*i*dh,  dh);
		}
		glEnd();
	}

	void drawCell(int i, int j){
		double dh = 1.0-offset;
		double wid = 2.0*dh/in_num;
		double hei = 2.0*dh/out_num;

		double r = 1.0/in_d_deg;
		for(double k=0.0; k<in_d_deg; k+=1){
			double x1 = (i + r*k)*wid - dh;
			double y1 =  j*hei - dh;
			double x2 = x1 + r*wid;
			double y2 = y1 + hei;
			double d = r*(k+0.5);
			MyColor c1 = (1.0 - d) * Color(i,j  ) + d * Color(i+1,j  );
			MyColor c2 = (1.0 - d) * Color(i,j+1) + d * Color(i+1,j+1);
			MyColor c = 0.5*(c1+c2);
			glColor3d(c.r, c.g, c.b);
			glRectd(x1, y1, x2, y2);
/*			
			//点対称なので,反対側にも表示
			x1 = (in_num - i - r*k)*wid - dh;
			y1 = (out_num -j)*hei - dh;
			x2 = x1 - r*wid;
			y2 = y1 - hei;
			glRectd(x1, y1, x2, y2);
*/		
		}
	}

	void drawSample(){
		int lambda = (700 - 380)/5.0; 
		double wid = (2.0-1.2*offset)/lambda;
		double hei = (2.0-1.2*offset)/out_num;
		double dh = 1.0-offset;
		for(int i=0; i<=lambda;i++){
			MyColor c1 = rgb->CIE_RGB(rgb->getXYZ(5*i+380,     0.2));
			MyColor c2 = rgb->CIE_RGB(rgb->getXYZ(5*(i+1)+380, 0.2));
			double r = 1.0/5.0;
			for(int k=0; k<5;k++){
				double x1 = (i+k*r)*wid - dh;
				double x2 = x1 + r;
				MyColor c = r*c2 + (1.0 - r)*c1;
				glColor3d(c.r, c.g, c.b);
				glRectd(x1, -1+dh, x2, 1.0);
			}
		}
	}

	void draw(){
		
		for(int i=0; i<in_num; i++)
			for(int j=0; j<out_num; j++)
				drawCell(i,j);
		draw_axis();
		
		//drawSample();
	}

	MyColor& Color(const int &i, const int &j){
		int a = (out_num+1)*i + j;
		return ColorMap[a];
	}

};
#endif