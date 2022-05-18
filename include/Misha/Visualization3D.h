/*
Copyright (c) 2021, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef VISUALIZATION_3D_INCLUDED
#define VISUALIZATION_3D_INCLUDED
#include "Misha/Visualization.h"
#include "Misha/Camera.h"
#include "Misha/Geometry.h"

#define NEW_CALL_BACK

namespace Misha
{
	template< typename DerivedViewableType >
	struct Viewable3D : public Viewable< DerivedViewableType >
	{
		using Viewable< DerivedViewableType >::screenWidth;
		using Viewable< DerivedViewableType >::screenHeight;
		using Viewable< DerivedViewableType >::callBacks;
		typedef struct Viewable< DerivedViewableType >::KeyboardCallBack KeyboardCallBack;

		enum CameraType
		{
			ORTHOGRAPHIC_CAMERA ,
			PERSPECTIVE_CAMERA
		};
		struct PerspectiveCameraInfo
		{
			Camera camera;
			double fovY , nearZ , farZ ;
			double zoom , scale;
		};
		struct OrthographicCameraInfo
		{
			Camera camera;
			double zoom , scale;
		};
		CameraType cameraType;
		PerspectiveCameraInfo perspectiveCameraInfo;
		OrthographicCameraInfo orthographicCameraInfo;
		Point3D< double > center;

		bool cameraForwardLight;
		GLfloat lightPosition[4];
		GLfloat lightAmbient[4] , lightDiffuse[4] , lightSpecular[4] , shapeSpecular[4] , shapeSpecularShininess;

		int oldX , oldY , newX , newY;
		bool rotating , scaling , panning;

		Viewable3D( void );
		void init( const Point3D< double > bBox[2] , unsigned int width , unsigned int height );

		void keyboardFunc( unsigned char key , int x , int y );
		void specialFunc( int key, int x, int y );
		void display( void );
		void mouseFunc( int button , int state , int x , int y );
		void motionFunc( int x , int y );

		bool select( int x , int y , Point3D< double > & p );

		virtual void idle( void ){}
		virtual void draw3D( void ) = 0;

		static void ToggleCameraType( DerivedViewableType *v , const char * )
		{
			if     ( v->cameraType==ORTHOGRAPHIC_CAMERA ) v->cameraType =  PERSPECTIVE_CAMERA;
			else if( v->cameraType== PERSPECTIVE_CAMERA ) v->cameraType = ORTHOGRAPHIC_CAMERA;
		}

		const Camera &camera( void ) const
		{
			if( cameraType==ORTHOGRAPHIC_CAMERA ) return orthographicCameraInfo.camera;
			else                                  return perspectiveCameraInfo.camera;
		}

		Camera &camera( void )
		{
			if( cameraType==ORTHOGRAPHIC_CAMERA ) return orthographicCameraInfo.camera;
			else                                  return perspectiveCameraInfo.camera;
		}

		float zoom( void ) const
		{
			if( cameraType==ORTHOGRAPHIC_CAMERA ) return orthographicCameraInfo.zoom;
			else                                  return perspectiveCameraInfo.zoom;
		}
	};


	////////////////
	// Viewable3D //
	////////////////

	template< typename DerivedViewableType >
	Viewable3D< DerivedViewableType >::Viewable3D( void )
	{
		cameraType = PERSPECTIVE_CAMERA;
		perspectiveCameraInfo.scale = orthographicCameraInfo.scale = 1.0;
		perspectiveCameraInfo.zoom = orthographicCameraInfo.zoom = 1.05;
		perspectiveCameraInfo.fovY = 30.;

		rotating = scaling = panning = false;

		cameraForwardLight = true;
		lightAmbient [0] = lightAmbient [1] = lightAmbient [2] = 0.25f , lightAmbient [3] = 1.f;
		lightDiffuse [0] = lightDiffuse [1] = lightDiffuse [2] = 0.70f , lightDiffuse [3] = 1.f;
		lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 1.00f , lightSpecular[3] = 1.f;
		shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 1.00f , shapeSpecular[3] = 1.f;
		shapeSpecularShininess = 128;

		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'C' , typename KeyboardCallBack::Modifiers() , "toggle camera" , ToggleCameraType ) );
	}

	template< typename DerivedViewableType >
	void Viewable3D< DerivedViewableType >::init( const Point3D< double > bBox[2] , unsigned int width , unsigned int height )
	{
		screenWidth = width , screenHeight = height;
		center = ( bBox[0] + bBox[1] ) / 2;
		double radius = Point3D< double >::Length( bBox[0]-center );
		perspectiveCameraInfo.scale = orthographicCameraInfo.scale = 1./radius;
		perspectiveCameraInfo.nearZ = 0.1*radius , perspectiveCameraInfo.farZ = 10.*radius;
		double tmp = ( sin( perspectiveCameraInfo.fovY * M_PI / 180. / 2. ) * width ) / height;
		if( tmp>1 ) THROW( "bad aspect ratio" );
		double fovX = ( asin( tmp ) * 2. * 180. / M_PI );

		// Given the unit cube, we would like to find the z-position such that the back face of the cube is completely visible:
		//    radius / (z-radius) = tan( fov/2 )
		// => z-radius = radius/tan( fov/2 )
		// => z = radius/tan( fov/2 ) + radius
		double zX = radius / tan( perspectiveCameraInfo.fovY * M_PI / 180. / 2. ) + radius;
		double zY = radius / tan( fovX * M_PI / 180. / 2. ) + radius;

		Point3D< double > position( center[0] , center[1] , center[2] + std::max< double >( zX , zY )+0.1 ) , forward( 0 , 0 , -1 ) , up( 0 , 1 , 0 );
		perspectiveCameraInfo.camera = Camera( position , forward , up );	
	}

	template< typename DerivedViewableType >
	void Viewable3D< DerivedViewableType >::display( void )
	{
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		if( cameraType==PERSPECTIVE_CAMERA )
		{
			GLint viewport[4];
			glGetIntegerv( GL_VIEWPORT , viewport );
			gluPerspective( perspectiveCameraInfo.fovY , (float)viewport[2] / viewport[3] , perspectiveCameraInfo.nearZ , perspectiveCameraInfo.farZ );
		}
		else if( cameraType==ORTHOGRAPHIC_CAMERA )
		{
			double ar = screenWidth/(double)screenHeight , ar_r = 1./ar;
			double zoom = orthographicCameraInfo.zoom;
			double scale = orthographicCameraInfo.scale;
			if( screenWidth>screenHeight ) glOrtho( -ar*zoom , ar*zoom , -zoom , zoom , -10./scale , 10./scale );
			else                           glOrtho( -zoom , zoom , -ar_r*zoom , ar_r*zoom , -10./scale , 10./scale );
		}

		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		if     ( cameraType== PERSPECTIVE_CAMERA )  perspectiveCameraInfo.camera.draw();
		else if( cameraType==ORTHOGRAPHIC_CAMERA ) orthographicCameraInfo.camera.draw();

		if( cameraForwardLight )
		{
			Point3D< float > d;
			if     ( cameraType== PERSPECTIVE_CAMERA ) d = Point3D< float >(  perspectiveCameraInfo.camera.up +  perspectiveCameraInfo.camera.right -  perspectiveCameraInfo.camera.forward*5 );
			else if( cameraType==ORTHOGRAPHIC_CAMERA ) d = Point3D< float >( orthographicCameraInfo.camera.up + orthographicCameraInfo.camera.right - orthographicCameraInfo.camera.forward*5 );
			lightPosition[0] = d[0] , lightPosition[1] = d[1] , lightPosition[2] = d[2];
			lightPosition[3] = 0;
		}
		glLightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER , GL_FALSE );
		glLightModeli( GL_LIGHT_MODEL_TWO_SIDE , GL_TRUE );
		glLightfv( GL_LIGHT0 , GL_AMBIENT , lightAmbient );
		glLightfv( GL_LIGHT0 , GL_DIFFUSE , lightDiffuse );
		glLightfv( GL_LIGHT0 , GL_SPECULAR , lightSpecular );
		glLightfv( GL_LIGHT0 , GL_POSITION , lightPosition );
		glEnable( GL_LIGHT0 );
		glEnable ( GL_LIGHTING );

		glColorMaterial( GL_FRONT_AND_BACK , GL_AMBIENT_AND_DIFFUSE );
		glEnable( GL_COLOR_MATERIAL );

		glEnable( GL_DEPTH_TEST );
		glMaterialfv( GL_FRONT_AND_BACK , GL_SPECULAR  , shapeSpecular );
		glMaterialf ( GL_FRONT_AND_BACK , GL_SHININESS , shapeSpecularShininess );
		glLightModeli( GL_LIGHT_MODEL_TWO_SIDE , GL_TRUE );

		static_cast< DerivedViewableType * >( this )->draw3D();
	}

	template< typename DerivedViewableType >
	void Viewable3D< DerivedViewableType >::mouseFunc( int button , int state , int x , int y )
	{
		newX = x ; newY = y;

		rotating = scaling = panning = false;
		if( button==GLUT_LEFT_BUTTON  )
		{
			if( state==GLUT_DOWN )
			{
				if( glutGetModifiers() & GLUT_ACTIVE_CTRL ) panning = true;
				else                                        rotating = true;
			}
		}
		else if( button==GLUT_RIGHT_BUTTON )
			scaling = true;
	}

	template< typename DerivedViewableType >
	void Viewable3D< DerivedViewableType >::motionFunc( int x , int y )
	{
		oldX = newX , oldY = newY , newX = x , newY = y;

		double zoom;
		if     ( cameraType== PERSPECTIVE_CAMERA ) zoom =  perspectiveCameraInfo.zoom;
		else if( cameraType==ORTHOGRAPHIC_CAMERA ) zoom = orthographicCameraInfo.zoom;

		int imageSize = std::min< int >( screenWidth , screenHeight );
		double rel_x = (newX - oldX) / (double)imageSize * 2;
		double rel_y = (newY - oldY) / (double)imageSize * 2;
		double pRight = -rel_x * zoom , pUp = rel_y * zoom;
		double sForward = rel_y*4;
		double rRight = rel_y , rUp = rel_x;

		if( cameraType==PERSPECTIVE_CAMERA )
		{
			Camera &camera = perspectiveCameraInfo.camera;
			if     ( rotating ) camera.rotateUp( rUp , center ) , camera.rotateRight( rRight , center );
			else if( scaling  ) zoom *= pow( 0.9 , (double)sForward );
			else if( panning  ) camera.translate( camera.right * pRight + camera.up * pUp );
		}
		else if( cameraType==ORTHOGRAPHIC_CAMERA )
		{
			Camera &camera = orthographicCameraInfo.camera;
			if     ( rotating ) camera.rotateUp( rUp , center ) , camera.rotateRight( rRight , center );
			else if( scaling  ) zoom *= pow( 0.9 , (double)sForward );
			else if( panning  ) camera.translate( camera.right * pRight + camera.up * pUp );
		}

		glutPostRedisplay();
	}

	template< typename DerivedViewableType >
	bool Viewable3D< DerivedViewableType >::select( int x , int y , Point3D< double > &p )
	{
		float z;
		glReadPixels( x , y , 1 , 1 , GL_DEPTH_COMPONENT , GL_FLOAT , &z );
		if( z<1 )
		{
			GLdouble mMatrix[16] , pMatrix[16];
			GLint viewport[4];
			glGetDoublev( GL_MODELVIEW_MATRIX , mMatrix );
			glGetDoublev( GL_PROJECTION_MATRIX , pMatrix );
			glGetIntegerv( GL_VIEWPORT , viewport );
			gluUnProject( (double)x , (double)y , (double)z , mMatrix , pMatrix , viewport , &p[0] , &p[1] , &p[2] );
			return true;
		}
		else return false;
	}

	template< typename DerivedViewableType >
	void Viewable3D< DerivedViewableType >::keyboardFunc( unsigned char key , int x , int y )
	{
		Camera &camera = cameraType==PERSPECTIVE_CAMERA ? perspectiveCameraInfo.camera : orthographicCameraInfo.camera;
		switch( key )
		{
		case 'q': camera.rotateUp( -M_PI/128 , center ) ; break;
		case 'Q': camera.rotateUp( -M_PI/  4 , center ) ; break;
		case 'w': camera.rotateUp(  M_PI/128 , center ) ; break;
		case 'W': camera.rotateUp(  M_PI/  4 , center ) ; break;
		case 'a': camera.rotateRight( -M_PI/128 , center ) ; break;
		case 'A': camera.rotateRight( -M_PI/  4 , center ) ; break;
		case 'z': camera.rotateRight(  M_PI/128 , center ) ; break;
		case 'Z': camera.rotateRight(  M_PI/  4 , center ) ; break;
		case 's': camera.rotateForward( -M_PI/128 , center ) ; break;
		case 'S': camera.rotateForward( -M_PI/  4 , center ) ; break;
		case 'x': camera.rotateForward(  M_PI/128 , center ) ; break;
		case 'X': camera.rotateForward(  M_PI/  4 , center ) ; break;
		}
	}

	template< typename DerivedViewableType >
	void Viewable3D< DerivedViewableType >::specialFunc( int key, int x, int y )
	{
//		double stepSize = 20. / ( screenWidth + screenHeight );
		double stepSize = (double)( screenWidth + screenHeight ) / 10000;
		if( glutGetModifiers()&GLUT_ACTIVE_CTRL ) stepSize /= 16;
		double panSize = stepSize*2;

		if( cameraType==PERSPECTIVE_CAMERA )
		{
			Camera &camera = perspectiveCameraInfo.camera;
			panSize /= perspectiveCameraInfo.scale;
			switch( key )
			{
			case Misha::KEY_UPARROW:    camera.translate(  camera.forward * panSize ) ; break;
			case Misha::KEY_DOWNARROW:  camera.translate( -camera.forward * panSize ) ; break;
			case Misha::KEY_LEFTARROW:  camera.translate(  camera.right   * panSize ) ; break;
			case Misha::KEY_RIGHTARROW: camera.translate( -camera.right   * panSize ) ; break;
			case Misha::KEY_PGUP:       camera.translate( -camera.up      * panSize ) ; break;
			case Misha::KEY_PGDN:       camera.translate(  camera.up      * panSize ) ; break;
			}
		}
		else if( cameraType==ORTHOGRAPHIC_CAMERA )
		{
			Camera &camera = orthographicCameraInfo.camera;
			switch( key )
			{
			case Misha::KEY_UPARROW:    orthographicCameraInfo.zoom *= 0.98f ; break;
			case Misha::KEY_DOWNARROW:  orthographicCameraInfo.zoom /= 0.98f ; break;
			case Misha::KEY_LEFTARROW:  camera.translate(  camera.right * panSize ) ; break;
			case Misha::KEY_RIGHTARROW: camera.translate( -camera.right * panSize ) ; break;
			case Misha::KEY_PGUP:       camera.translate( -camera.up    * panSize ) ; break;
			case Misha::KEY_PGDN:       camera.translate(  camera.up    * panSize ) ; break;
			}
		}
		glutPostRedisplay();
	}
}
#endif // VISUALIZATION_3D_INCLUDED