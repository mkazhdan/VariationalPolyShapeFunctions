/*
Copyright (c) 2018, Michael Kazhdan
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

#ifndef VISUALIZATION_INCLUDED
#define VISUALIZATION_INCLUDED
#include <algorithm>
#include <sstream>
#include <GL/glew.h>
#include <GL/glut.h>
#include <vector>
#include <sys/timeb.h>
#include <functional>
#include <algorithm>
#include "Misha/Image.h"
#include "Misha/Array.h"
#include "Misha/Exceptions.h"
#include "Misha/Miscellany.h"

#define NEW_CALL_BACK

namespace Misha
{
	static const int KEY_UPARROW    = 101;
	static const int KEY_DOWNARROW	= 103;
	static const int KEY_LEFTARROW	= 100;
	static const int KEY_RIGHTARROW	= 102;
	static const int KEY_PGUP		= 104;
	static const int KEY_PGDN		= 105;
	static const int KEY_CTRL_C     =   3;
	static const int KEY_BACK_SPACE =   8;
	static const int KEY_ENTER      =  13;
	static const int KEY_ESC        =  27;

	struct Font
	{
		std::string fontName;
		void *font;
		unsigned int fontHeight;

		Font( std::string fn , void *f , unsigned int fh ) : fontName(fn) , font(f) , fontHeight(fh) {}

		static const Font Fonts[];
		static unsigned int FontNum();
	};

	const Font Font::Fonts[] =
	{
		Font( "Bitmap 8x13" , GLUT_BITMAP_8_BY_13 , 13 ) ,
		Font( "Bitmap 9x15" , GLUT_BITMAP_9_BY_15 , 15 ) ,
		Font( "Helvetica 10" , GLUT_BITMAP_HELVETICA_10 , 10 ) ,
		Font( "Helvetica 12" , GLUT_BITMAP_HELVETICA_12 , 12 ) ,
		Font( "Helvetica 18" , GLUT_BITMAP_HELVETICA_18 , 18 ) ,
		Font( "Times-Roman 10" , GLUT_BITMAP_TIMES_ROMAN_10 , 10 ) ,
		Font( "Times-Roman 24" , GLUT_BITMAP_TIMES_ROMAN_24 , 24 ) ,
	};
	unsigned int Font::FontNum( void ){ return sizeof( Fonts ) / sizeof( Font ); }

	template< typename DerivedViewableType >
	struct Viewable
	{
		struct KeyboardCallBack
		{
#ifdef NEW_CALL_BACK
			struct Modifiers
			{
				bool alt , ctrl;
				Modifiers( bool a=false , bool c=false ) : alt(a) , ctrl(c) {};
				bool operator == ( const Modifiers &m ) const { return alt==m.alt && ctrl==m.ctrl; }
				bool operator != ( const Modifiers &m ) const { return alt!=m.alt || ctrl!=m.ctrl; }
			};
#endif // NEW_CALL_BACK
			char key;
			char prompt[1024];
			char description[1024];
			void (*callBackFunction)( DerivedViewableType* , const char* );
			DerivedViewableType* viewable;
#ifdef NEW_CALL_BACK
			Modifiers modifiers;
			KeyboardCallBack( DerivedViewableType* viewable , char key , Modifiers modifiers , const char* description ,                      void ( *callBackFunction )( DerivedViewableType* , const char* ) );
			KeyboardCallBack( DerivedViewableType* viewable , char key , Modifiers modifiers , const char* description , const char* prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) );
#else // !NEW_CALL_BACK
			KeyboardCallBack( DerivedViewableType* viewable , char key , const char* description ,                      void ( *callBackFunction )( DerivedViewableType* , const char* ) );
			KeyboardCallBack( DerivedViewableType* viewable , char key , const char* description , const char* prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) );
#endif // NEW_CALL_BACK
		};
		struct Viewer
		{
			static DerivedViewableType *viewable;
			static void Run( DerivedViewableType* viewable , int argc , char* argv[] , const char* windowName="" );
			static void Idle             ( void );
			static void KeyboardFunc     ( unsigned char key , int x , int y );
			static void KeyboardUpFunc   ( unsigned char key , int x , int y );
			static void SpecialFunc      ( int key, int x, int y );
			static void SpecialUpFunc    ( int key, int x, int y );
			static void Display          ( void );
			static void Reshape          ( int w , int h );
			static void MouseFunc        ( int button , int state , int x , int y );
			static void MotionFunc       ( int x , int y );
			static void PassiveMotionFunc( int x , int y );
		};

	protected:
		const double _MIN_FPS_TIME = 0.5;
		Miscellany::Timer _fpsTimer;
//		double _lastFPSTime;
		int _lastFPSCount;
		double _fps;
		int _currentFrame , _totalFrames;
		bool _exitAfterSnapshot , _exitAfterVideo;
	public:
		int screenWidth , screenHeight;
		void *font , *promptFont;
		int fontHeight , promptFontHeight;
		bool showHelp , showInfo , showFPS;
		void (*promptCallBack)( DerivedViewableType* , const char* );
		char promptString[1024];
		int promptLength;
		char* snapshotName;
		char* videoHeader;
		bool flushImage;
		std::function< void ( void ) > quitFunction;

		std::vector< KeyboardCallBack > callBacks;
		std::vector< char* > info;
		Viewable( void );
		void setFont( unsigned int idx );
		void setSnapshot( const char* sName , bool exitAfterSnapshot=true );
		void setVideo( const char* sName , int frames , bool exitAfterVideo=true );
		virtual void display( void ) {}
		virtual void idle( void ) {}
		virtual void keyboardFunc( unsigned char key , int x , int y ) {}
		virtual void keyboardUpFunc( unsigned char key , int x , int y ) {}
		virtual void specialFunc( int key, int x, int y ) {}
		virtual void specialUpFunc( int key, int x, int y ) {}
		virtual void mouseFunc( int button , int state , int x , int y ) {}
		virtual void motionFunc( int x , int y ) {}
		virtual void passiveMotionFunc( int x , int y ) {}
		virtual void reshape( int w , int h )
		{
			screenWidth = w , screenHeight = h;
			glViewport( 0 , 0 , screenWidth , screenHeight );
		}

		void Idle             ( void );
		void KeyboardFunc     ( unsigned char key , int x , int y );
		void KeyboardUpFunc   ( unsigned char key , int x , int y );
		void SpecialFunc      ( int key, int x, int y );
		void SpecialUpFunc    ( int key, int x, int y );
		void Display          ( void );
		void Reshape          ( int w , int h );
		void MouseFunc        ( int button , int state , int x , int y );
		void MotionFunc       ( int x , int y );
		void PassiveMotionFunc( int x , int y );

#if 0
#pragma message( "[WARNING] Disabling exit" )
		static void           ExitCallBack( DerivedViewableType*   , const char* ){                  ; }
		static void           QuitCallBack( DerivedViewableType* v , const char* ){ v->quitFunction(); }
#else
		static void           ExitCallBack( DerivedViewableType*   , const char* ){                     exit( 0 ); }
		static void           QuitCallBack( DerivedViewableType* v , const char* ){ v->quitFunction() , exit( 0 ); }
#endif
		static void      ToggleFPSCallBack( DerivedViewableType* v , const char* ){ v->showFPS  = !v->showFPS ; }
		static void     ToggleHelpCallBack( DerivedViewableType* v , const char* ){ v->showHelp = !v->showHelp; }
		static void     ToggleInfoCallBack( DerivedViewableType* v , const char* ){ v->showInfo = !v->showInfo; }
		static void SetFrameBufferCallBack( DerivedViewableType* v , const char* prompt );

		static void WriteLeftString( int x , int y , void* font , const char* format , ... );
		static int StringWidth( void* font , const char* format , ... );
		void writeLeftString( int x , int y , const char* format , ... ) const;
		void writeRightString( int x , int y , const char* format , ... ) const;
		void writeCenterString( int x , int y , const char* format , ... ) const;
		int stringWidth( const char* format , ... ) const;

		void saveFrameBuffer( const char* fileName , int whichBuffer=GL_BACK );

		void setPromptCallBack( const char *prompt , void (*callBackFunction)( DerivedViewableType* , const char* ) );
#if 0
		if( strlen( callBacks[i].prompt ) )
		{
			sprintf( promptString , "%s: " , callBacks[i].prompt );
			promptLength = int( strlen( promptString ) );
			promptCallBack = callBacks[i].callBackFunction;
		}
		else (*callBacks[i].callBackFunction)( (DerivedViewableType*)this , NULL );
#endif


#ifdef NEW_CALL_BACK
		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char *description ,                      void ( *callBackFunction )( DerivedViewableType* , const char* ) );
		void addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char *description , const char *prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) );
#endif // NEW_CALL_BACK
		void addCallBack( char key , const char *description ,                      void ( *callBackFunction )( DerivedViewableType* , const char* ) );
		void addCallBack( char key , const char *description , const char *prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) );
	};

	//////////////////////
	// KeyboardCallBack //
	//////////////////////
	template< typename DerivedViewableType >
#ifdef NEW_CALL_BACK
	Viewable< DerivedViewableType >::KeyboardCallBack::KeyboardCallBack( DerivedViewableType* viewable , char key , Modifiers modifiers , const char* description , void (*callBackFunction)( DerivedViewableType* , const char* ) )
#else // !NEW_CALL_BACK
	Viewable< DerivedViewableType >::KeyboardCallBack::KeyboardCallBack( DerivedViewableType* viewable , char key , const char* description , void (*callBackFunction)( DerivedViewableType* , const char* ) )
#endif // NEW_CALL_BACK
	{
#ifdef NEW_CALL_BACK
		this->modifiers = modifiers;
#endif // NEW_CALL_BACK
		this->viewable = viewable;
		this->key = key;
		strcpy( this->description , description );
		prompt[0] = 0;
		this->callBackFunction = callBackFunction;
	}

	template< typename DerivedViewableType >
#ifdef NEW_CALL_BACK
	Viewable< DerivedViewableType >::KeyboardCallBack::KeyboardCallBack( DerivedViewableType* viewable , char key , Modifiers modifiers , const char* description , const char* prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) )
#else // !NEW_CALL_BACK
	Viewable< DerivedViewableType >::KeyboardCallBack::KeyboardCallBack( DerivedViewableType* viewable , char key , const char* description , const char* prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) )
#endif // NEW_CALL_BACK
	{
#ifdef NEW_CALL_BACK
		this->modifiers = modifiers;
#endif // NEW_CALL_BACK
		this->viewable = viewable;
		this->key = key;
		strcpy( this->description , description );
		strcpy( this->prompt , prompt );
		this->callBackFunction = callBackFunction;
	}

	//////////////////////
	// Viewable::Viewer //
	//////////////////////

	template< typename DerivedViewableType > DerivedViewableType* Viewable< DerivedViewableType >::Viewer::viewable = NULL;

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::Viewer::Run( DerivedViewableType* v , int argc , char* argv[] , const char* windowName )
	{
		viewable = v;
		glutInitDisplayMode( GLUT_DEPTH | GLUT_RGB | GLUT_DOUBLE );
		glutInitWindowSize( viewable->screenWidth , viewable->screenHeight );
		glutInit( &argc , argv );
		glutCreateWindow( windowName );

		if( glewInit()!=GLEW_OK ) ERROR_OUT( "glewInit failed" );
		glutIdleFunc         ( Idle );
		glutDisplayFunc      ( Display );
		glutReshapeFunc      ( Reshape );
		glutMouseFunc        ( MouseFunc );
		glutMotionFunc       ( MotionFunc );
		glutPassiveMotionFunc( PassiveMotionFunc );
		glutKeyboardFunc     ( KeyboardFunc );
		glutKeyboardUpFunc   ( KeyboardUpFunc );
		glutSpecialFunc      ( SpecialFunc );
		glutSpecialUpFunc    ( SpecialUpFunc );

		glutMainLoop();
	}

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::Idle( void ){ viewable->Idle(); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::KeyboardFunc( unsigned char key , int x , int y ){ viewable->KeyboardFunc( key , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::KeyboardUpFunc( unsigned char key , int x , int y ){ viewable->KeyboardUpFunc( key , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::SpecialFunc( int key , int x , int y ){ viewable->SpecialFunc( key , x ,  y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::SpecialUpFunc( int key , int x , int y ){ viewable->SpecialUpFunc( key , x ,  y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::Display( void ){ viewable->Display(); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::Reshape( int w , int h ){ viewable->Reshape( w , h ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::MouseFunc( int button , int state , int x , int y ){ viewable->MouseFunc( button , state , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::MotionFunc( int x , int y ){ viewable->MotionFunc( x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Viewer::PassiveMotionFunc( int x , int y ){ viewable->PassiveMotionFunc( x , y ); }

	//////////////
	// Viewable //
	//////////////
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::Reshape( int w , int h ){ reshape( w , h ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::MouseFunc( int button , int state , int x , int y ){ mouseFunc( button , state , x , y ); }
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::MotionFunc( int x , int y ){ motionFunc( x , y );}
	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::PassiveMotionFunc( int x , int y ){ passiveMotionFunc( x , y );}
	template< typename DerivedViewableType > 
	void Viewable< DerivedViewableType >::Idle( void )
	{
		if( snapshotName )
		{
			if( flushImage )
			{
				flushImage = false;
				glutPostRedisplay();
				return;
			}
			else
			{
				saveFrameBuffer( snapshotName , GL_FRONT );
				delete[] snapshotName;
				snapshotName = NULL;
				if( _exitAfterSnapshot ) exit( 0 );
			}
		}
		else if( videoHeader && _currentFrame<_totalFrames )
		{
			char snapshotName[512];
//			sprintf( snapshotName , "%s.%04d.jpg" , videoHeader , _currentFrame );
			sprintf( snapshotName , "%s.%04d.png" , videoHeader , _currentFrame );
			saveFrameBuffer( snapshotName , GL_FRONT );
			_currentFrame++;
			if( _currentFrame==_totalFrames && _exitAfterVideo ) exit( 0 );
		}
		idle();
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setPromptCallBack( const char *prompt , void (*callBackFunction)( DerivedViewableType* , const char* ) )
	{
		sprintf( promptString , "%s: " , prompt );
		promptLength = (int)strlen( promptString );
		promptCallBack = callBackFunction;
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::KeyboardFunc( unsigned char key , int x , int y )
	{
		if( promptCallBack )
		{
			size_t len = strlen( promptString );
			if( key==KEY_BACK_SPACE )
			{
				if( len>promptLength ) promptString[len-1] = 0;
			}
			else if( key==KEY_ENTER )
			{
				promptCallBack( (DerivedViewableType*)this , promptString+promptLength );
				promptString[0] = 0;
				promptLength = 0;
				promptCallBack = NULL;
			}
			else if( key==KEY_CTRL_C )
			{
				promptString[0] = 0;
				promptLength = 0;
				promptCallBack = NULL;
			}
			else if( key>=32 && key<=126 ) // ' ' to '~'
			{
				promptString[ len ] = key;
				promptString[ len+1 ] = 0;
			}
			glutPostRedisplay();
			return;
		}
		switch( key )
		{
		case KEY_CTRL_C:
			exit( 0 );
			break;
		default:
#ifdef NEW_CALL_BACK
		{
			int m = glutGetModifiers();
			typename KeyboardCallBack::Modifiers modifiers( m & GLUT_ACTIVE_ALT , m & GLUT_ACTIVE_CTRL );
			for( int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].key==key && callBacks[i].modifiers==modifiers )
			{
				if( strlen( callBacks[i].prompt ) ) setPromptCallBack( callBacks[i].prompt , callBacks[i].callBackFunction );
				else (*callBacks[i].callBackFunction)( (DerivedViewableType*)this , NULL );
				break;
			}
		}
#else // !NEW_CALL_BACK
			for( int i=0 ; i<callBacks.size() ; i++ ) if( callBacks[i].key==key )
			{
				if( strlen( callBacks[i].prompt ) ) setPromptCallBack( callBacks[i].prompt , callBacks[i].callBackFunction );
				else (*callBacks[i].callBackFunction)( (DerivedViewableType*)this , NULL );
				break;
			}
#endif // NEW_CALL_BACK
		}
		keyboardFunc( key , x , y );
		glutPostRedisplay();
	}

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::KeyboardUpFunc( unsigned char key , int x , int y ){ keyboardUpFunc( key , x , y ); }

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::SpecialFunc( int key , int x , int y ){ specialFunc( key , x , y ); }

	template< typename DerivedViewableType > void Viewable< DerivedViewableType >::SpecialUpFunc( int key , int x , int y ){ specialUpFunc( key , x , y ); }

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::Display( void )
	{
		glClearColor( 1 , 1 , 1 , 1 );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

		display();

		_lastFPSCount++;
#if 1
		double t = _fpsTimer.elapsed();
		if( t > _MIN_FPS_TIME )
		{
			_fps = (double)_lastFPSCount / t;
			_lastFPSCount = 0;
			_fpsTimer.reset();
		}
#else
		double t = Time();
		if( t-_lastFPSTime > _MIN_FPS_TIME )
		{
			_fps = (double)_lastFPSCount / (t-_lastFPSTime);
			_lastFPSCount = 0;
			_lastFPSTime = t;
		}
#endif
		if( showFPS ) writeRightString( 5 , screenHeight - fontHeight - 5 , "%d x %d @ %.2f" , screenWidth , screenHeight , _fps );

		GLboolean writeMask;
		glGetBooleanv( GL_DEPTH_WRITEMASK , &writeMask );
		glDepthMask( GL_FALSE );

		glDisable( GL_LIGHTING );
		int offset = fontHeight/2;
		if( showHelp )
		{
#ifdef NEW_CALL_BACK
			auto CallBackString = [&]( int i )
			{
				std::stringstream stream;
				stream << "\'" << callBacks[i].key << "\'";
				if( callBacks[i].modifiers.ctrl ) stream << "+[CTRL]";
				if( callBacks[i].modifiers.alt ) stream << "+[ALT]";
				stream << ": " << callBacks[i].description;
				return stream.str();
			};
#endif // NEW_CALL_BACK
			{
				GLint vp[4];
				glGetIntegerv( GL_VIEWPORT , vp );

				glMatrixMode( GL_PROJECTION );
				glPushMatrix();
				glLoadIdentity();
				glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

				glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
				glLoadIdentity();

				int x=0 , y = offset;
#ifdef NEW_CALL_BACK
				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) ) x = std::max< int >( x , stringWidth( CallBackString(i).c_str() ) ) , y += fontHeight + offset;
#else // !NEW_CALL_BACK
				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) ) x = std::max< int >( x , stringWidth( "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) ) , y += fontHeight + offset;
#endif // NEW_CALL_BACK

				GLint srcAlpha , dstAlpha;
				glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
				glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );

				glEnable( GL_BLEND );
				glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
				glBegin( GL_QUADS );
				glColor4f( 1.f , 1.f , 1.f , 0.5f );
				glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
				glEnd();
				glBlendFunc( srcAlpha , dstAlpha );
				glDisable( GL_BLEND );
				glDisable( GL_DEPTH_TEST );
				glLineWidth( 2.f );
				glBegin( GL_LINE_LOOP );
				glColor4f( 0.f , 0.f , 0.f , 1.f );
				glVertex2i( screenWidth-5 , 0 ) , glVertex2i( screenWidth-(x+15) , 0 ) , glVertex2i( screenWidth-(x+15) , y ) , glVertex2i( screenWidth-5 , y );
				glEnd();

				glMatrixMode( GL_PROJECTION );
				glPopMatrix();

				glMatrixMode( GL_MODELVIEW );
				glPopMatrix();
			}

			{
				int y = offset , width = 0;

				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
#ifdef NEW_CALL_BACK
					width = std::max< int >( width , stringWidth( CallBackString(i).c_str() ) );
#else // !NEW_CALL_BACK
					width = std::max< int >( width , stringWidth( "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) );
#endif // NEW_CALL_BACK
				for( int i=0 ; i<callBacks.size() ; i++ ) if( strlen( callBacks[i].description ) )
#ifdef NEW_CALL_BACK
					writeLeftString( screenWidth - 10 - width , y , CallBackString(i).c_str() ) , y += fontHeight + offset;
#else // !NEW_CALL_BACK
					writeLeftString( screenWidth - 10 - width , y , "\'%c\': %s" , callBacks[i].key , callBacks[i].description ) , y += fontHeight + offset;
#endif // NEW_CALL_BACK
			}
		}
		if( showInfo && info.size() )
		{
			{
				GLint vp[4];
				glGetIntegerv( GL_VIEWPORT , vp );

				glMatrixMode( GL_MODELVIEW );
				glPushMatrix();
				glLoadIdentity();
				glMatrixMode( GL_PROJECTION );
				glPushMatrix();
				glLoadIdentity();
				glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );
				int x=0 , y = offset;
				for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) ) x = std::max< int >( x , glutBitmapLength( font , (unsigned char*) info[i] ) ) , y += fontHeight + offset;
				glEnable( GL_BLEND );
				GLint srcAlpha , dstAlpha;
				glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
				glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
				glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
				glBegin( GL_QUADS );
				glColor4f( 1.f , 1.f , 1.f , 0.5f );
				glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
				glEnd();
				glBlendFunc( srcAlpha , dstAlpha );
				glDisable( GL_BLEND );
				glDisable( GL_DEPTH_TEST );
				glLineWidth( 2.f );
				glBegin( GL_LINE_LOOP );
				glColor4f( 0.f , 0.f , 0.f , 1.f );
				glVertex2i( 5 , 0 ) , glVertex2i( x+15 , 0 ) , glVertex2i( x+15 , y ) , glVertex2i( 5 , y );
				glEnd();

				glMatrixMode( GL_PROJECTION );
				glPopMatrix();

				glMatrixMode( GL_MODELVIEW );
				glPopMatrix();
			}
			{
				int y = offset;
				for( int i=0 ; i<info.size() ; i++ ) if( strlen( info[i] ) )
					writeLeftString( 10 , y , "%s" , info[i] ) , y += fontHeight + offset;
			}
		}
		if( strlen( promptString ) )
		{
			void* _font = font;
			int _fontHeight = fontHeight;
			font = promptFont;
			fontHeight = promptFontHeight;

			int sw = StringWidth ( font , promptString );
			glColor4f( 1.f , 1.f , 1.f , 0.5 );
			glEnable( GL_BLEND );
			GLint srcAlpha , dstAlpha;
			glGetIntegerv( GL_BLEND_SRC_ALPHA , &srcAlpha );
			glGetIntegerv( GL_BLEND_DST_ALPHA , &dstAlpha );
			glBlendFunc( GL_SRC_ALPHA , GL_ONE_MINUS_SRC_ALPHA );
			glBegin( GL_QUADS );
			{
				glVertex2i(     0 , screenHeight              );
				glVertex2i( sw+20 , screenHeight              );
				glVertex2i( sw+20 , screenHeight-fontHeight*2 );
				glVertex2i(     0 , screenHeight-fontHeight*2 );
			}
			glEnd();
			glBlendFunc( srcAlpha , dstAlpha );
			glDisable( GL_BLEND );
			glColor4f( 0.f , 0.f , 0.f , 1.f );
			glLineWidth( 2.f );
			glBegin( GL_LINE_LOOP );
			{
				glVertex2i(     0 , screenHeight              );
				glVertex2i( sw+20 , screenHeight              );
				glVertex2i( sw+20 , screenHeight-fontHeight*2 );
				glVertex2i(     0 , screenHeight-fontHeight*2 );
			}
			glEnd();
			writeLeftString( 10 , screenHeight-fontHeight-fontHeight/2 , promptString );
			font = _font;
			fontHeight = _fontHeight;
		}
		if( writeMask ) glDepthMask( GL_TRUE );
		glutSwapBuffers();
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::WriteLeftString( int x , int y , void* font , const char* format , ... )
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}

		GLint vp[4];
		glGetIntegerv( GL_VIEWPORT , vp );

		glMatrixMode( GL_PROJECTION );
		glPushMatrix();
		glLoadIdentity();
		glOrtho( vp[0] , vp[2] , vp[1] , vp[3] , 0 , 1 );

		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();
		glLoadIdentity();

		GLint matrixMode;
		glGetIntegerv( GL_MATRIX_MODE , &matrixMode );
		int depth = glIsEnabled( GL_DEPTH_TEST );
		int lighting = glIsEnabled( GL_LIGHTING );
		glDisable( GL_DEPTH_TEST );
		glDisable( GL_LIGHTING );
		glColor4f( 0 , 0 , 0 , 1 );
		glRasterPos2i( x , y );
		int len = int( strlen( str ) );
		for( int i=0 ; i<len ; i++ ) glutBitmapCharacter( font , str[i] );
		if( depth ) glEnable( GL_DEPTH_TEST );
		if( lighting ) glEnable( GL_LIGHTING );

		glMatrixMode( GL_PROJECTION );
		glPopMatrix();

		glMatrixMode( GL_MODELVIEW );
		glPopMatrix();

		glMatrixMode( matrixMode );
	}

	template< typename DerivedViewableType >
	int Viewable< DerivedViewableType >::StringWidth( void* font , const char* format , ... )
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		return glutBitmapLength( font , (unsigned char*) str );
	}

	template< typename DerivedViewableType >
	int Viewable< DerivedViewableType >::stringWidth( const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		return glutBitmapLength( font , (unsigned char*) str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeLeftString( int x , int y , const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( x , y , font , str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeRightString( int x , int y , const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( screenWidth-x-glutBitmapLength( font , (unsigned char*) str ) , y , font  ,str );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::writeCenterString( int x , int y , const char* format , ... ) const
	{
		static char str[1024];
		{
			va_list args;
			va_start( args , format );
			vsprintf( str , format , args );
			va_end( args );
		}
		WriteLeftString( x-glutBitmapLength( font , (unsigned char*) str )/2 , y , font  ,str );
	}


	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::saveFrameBuffer( const char* fileName , int whichBuffer )
	{
		Pointer( float ) pixels = AllocPointer< float >( sizeof(float) * 3 * screenWidth * screenHeight );
		Pointer( unsigned char ) _pixels = AllocPointer< unsigned char >( sizeof(unsigned char) * 3 * screenWidth * screenHeight );
		glReadBuffer( whichBuffer );
		glReadPixels( 0 , 0 , screenWidth , screenHeight , GL_RGB , GL_FLOAT , pixels );
		for( int j=0 ; j<screenHeight ; j++ ) for( int i=0 ; i<screenWidth ; i++ ) for( int c=0 ; c<3 ; c++ )
		{
			int ii = int( pixels[ c + i * 3 + ( screenHeight - 1 - j ) * screenWidth * 3 ]*256 );
			if( ii<  0 ) ii =   0;
			if( ii>255 ) ii = 255;
			_pixels[ c + i * 3 + j * screenWidth * 3 ] = (unsigned char)ii;
		}
		FreePointer( pixels );
		ImageWriter::Write( fileName , _pixels , screenWidth , screenHeight , 3 , 95 );
		FreePointer( _pixels );
	}

#ifdef NEW_CALL_BACK
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char *description , void ( *callBackFunction )( DerivedViewableType* , const char* ) )
	{
		callBacks.push_back( KeyboardCallBack( this , key , modifiers , description , callBackFunction ) );
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , typename KeyboardCallBack::Modifiers modifiers , const char *description , const char *prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) )
	{
		callBacks.push_back( KeyboardCallBack( this , key , modifiers , description , prompt , callBackFunction ) );
	}
#endif // NEW_CALL_BACK
	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , const char *description , void ( *callBackFunction )( DerivedViewableType* , const char* ) )
	{
#ifdef NEW_CALL_BACK
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , key , KeyboardCallBack::Modifiers() , description , callBackFunction ) );
#else // !NEW_CALL_BACK
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , key , description , callBackFunction ) );
#endif // NEW_CALL_BACK
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::addCallBack( char key , const char *description , const char *prompt , void ( *callBackFunction )( DerivedViewableType* , const char* ) )
	{
#ifdef NEW_CALL_BACK
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , key , KeyboardCallBack::Modifiers() , description , prompt , callBackFunction ) );
#else // !NEW_CALL_BACK
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , key , description , prompt , callBackFunction ) );
#endif // NEW_CALL_BACK
	}

	template< typename DerivedViewableType >
	Viewable< DerivedViewableType >::Viewable( void )
	{
		quitFunction = []( void ){};
#ifdef NEW_CALL_BACK
		typename KeyboardCallBack::Modifiers modifiers;
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , KEY_ESC    , modifiers , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , KEY_CTRL_C , modifiers , "" , ExitCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'F' , modifiers , "toggle fps"  , ToggleFPSCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'H' , modifiers , "toggle help" , ToggleHelpCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'I' , modifiers , "toggle info" , ToggleInfoCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'i' , modifiers , "save frame buffer" , "Ouput image" , SetFrameBufferCallBack ) );
#else // !NEW_CALL_BACK
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , KEY_ESC    , "" , QuitCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , KEY_CTRL_C , "" , ExitCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'F' , "toggle fps"  , ToggleFPSCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'H' , "toggle help" , ToggleHelpCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'I' , "toggle info" , ToggleInfoCallBack ) );
		callBacks.push_back( KeyboardCallBack( (DerivedViewableType*)this , 'i' , "save frame buffer" , "Ouput image" , SetFrameBufferCallBack ) );
#endif // NEW_CALL_BACK
		snapshotName = videoHeader = NULL;
		flushImage = false;
		showHelp = showInfo = showFPS = true;
		_exitAfterSnapshot = _exitAfterVideo = false;
		screenWidth = screenHeight = 512;
		font = GLUT_BITMAP_HELVETICA_12;
		fontHeight = 12;
		promptFont = GLUT_BITMAP_TIMES_ROMAN_24;
		promptFontHeight = 24;
		promptCallBack = NULL;
		strcpy( promptString , "" );
		promptLength = 0;

//		_lastFPSTime = Time();
		_lastFPSCount = 0;
		_fps = 0;
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setFont( unsigned int idx )
	{
		if( idx>=Font::FontNum() ) WARN( "Font index out of bounds: " , idx , " < " , Font::FontNum() );
		else font = Font::Fonts[idx].font , fontHeight = Font::Fonts[idx].fontHeight;
	}


	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setSnapshot( const char* sName , bool exitAfterSnapshot )
	{
		_exitAfterSnapshot = exitAfterSnapshot;
		snapshotName = new char[ strlen( sName ) + 1 ];
		strcpy( snapshotName , sName );
		showHelp = showInfo = showFPS = false;
		flushImage = true;
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::setVideo( const char* vHeader , int frames , bool exitAfterVideo )
	{
		_exitAfterVideo = exitAfterVideo;
		videoHeader = new char[ strlen( vHeader ) + 1 ];
		strcpy( videoHeader , vHeader );
		_currentFrame = 0;
		_totalFrames = frames;
	}

	template< typename DerivedViewableType >
	void Viewable< DerivedViewableType >::SetFrameBufferCallBack( DerivedViewableType* v , const char* prompt )
	{
		if( prompt )
		{
			v->snapshotName = new char[ strlen(prompt)+1 ];
			strcpy( v->snapshotName , prompt );
			v->flushImage = true;
		}
	}
}
#endif // VISUALIZATION_INCLUDED