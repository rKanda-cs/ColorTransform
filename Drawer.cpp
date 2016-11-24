#include"Drawer.h"
#include<glut.h>

/* 文字描画*/
void drawBitmapString(void *font, string s, double x, double y)
{
	glPushAttrib(GL_CURRENT_BIT);
	glRasterPos2d(x, y);		/* 描画位置 */
	/* ビットマップ文字列の描画 */
	for(unsigned int i=0; i<s.size() ;i++)
		glutBitmapCharacter(font, s[i]);

  glPopAttrib();
}

MyColor operator+ (MyColor &a, MyColor &b){
	return MyColor(a.r+b.r, a.g+b.g, a.b+b.b);
};

MyColor operator- (MyColor &a, MyColor &b){
	return MyColor(a.r-b.r, a.g-b.g, a.b-b.b);
}; 


MyColor operator* (const double &a, const MyColor &b){
	return MyColor(a*b.r, a*b.g, a*b.b);
}


MyColor operator/ (const MyColor &b, const double &a){
	return MyColor(b.r/a, b.g/a, b.b/a);
}

bool operator== (MyColor &a, MyColor &b){
	return (a.r=b.r) && (a.g == b.g) && (a.b == b.b);
}


namespace ColorList{
	MyColor getColor(int i){
		i = i % 25 + 1;
		double red   = 0.5 * (int)(i/9);
		double green = 0.5 * ((i/9)%3);
		double blue  = 0.5 * (i%3);
		return MyColor( red, green ,blue); 
	}
};

unsigned int DataSet<double>::ID = 0;


