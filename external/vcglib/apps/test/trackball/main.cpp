#include <iostream>
using namespace std;

#ifdef _WIN32
#include<windows.h>
#endif

#include <SDL/SDL.h>

//#include <GL/glew.h>
#include <wrap/gui/trackball.h>
#include <GL/glut.h>

using namespace vcg;

bool fullscreen = false;
//int width =1024;
//int height = 768;
int width = 800;
int height = 600;



SDL_Surface *screen = NULL;

bool init() {
  
  if(SDL_Init(SDL_INIT_VIDEO) != 0) {
    return false;
  }

  const SDL_VideoInfo *info = SDL_GetVideoInfo();
  int bpp = info->vfmt->BitsPerPixel;

  SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 5);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

  int flags = SDL_OPENGL;
  if(fullscreen) 
    flags |= SDL_FULLSCREEN;

  screen = SDL_SetVideoMode(width, height, bpp, flags);
  if(!screen) {
    return false;
  }
  
  SDL_WM_SetIcon(SDL_LoadBMP("inspector.bmp"), NULL);
  SDL_WM_SetCaption(" Inspector", "Inspector");


  glDisable(GL_DITHER);
  glShadeModel(GL_SMOOTH);
  glHint( GL_FOG_HINT, GL_NICEST );
  glEnable(GL_DEPTH_TEST);
  glDepthFunc( GL_LEQUAL );
  glDisable(GL_LIGHTING); 

  return true;
}




int main(int argc, unsigned short **argv) {  
  if(!init()) return -1;
  glewInit();
  
  Trackball trackball;

  int quit = 0;
  int x, y;
  SDL_Event event;  
  while( !quit ) {                
    while( SDL_PollEvent( &event ) ){                        
      switch( event.type ) {
        case SDL_QUIT:  quit = 1; break; 

        case SDL_KEYDOWN:   
          switch(event.key.keysym.sym) {
            case SDLK_RSHIFT:
            case SDLK_LSHIFT:
              trackball.ButtonDown(Trackball::KEY_SHIFT);
              break;

            case SDLK_RCTRL:
            case SDLK_LCTRL:
              trackball.ButtonDown(Trackball::KEY_CTRL);
              break;

            case SDLK_RALT:
            case SDLK_LALT:
              trackball.ButtonDown(Trackball::KEY_ALT);
              break;
          }
          break;

        case SDL_KEYUP:
          switch(event.key.keysym.sym) {
            case SDLK_q: exit(0); break;            
            
            case SDLK_RSHIFT:
            case SDLK_LSHIFT:
              trackball.ButtonUp(Trackball::KEY_SHIFT);
              break;

            case SDLK_RCTRL:
            case SDLK_LCTRL:
              trackball.ButtonUp(Trackball::KEY_CTRL);
              break;

            case SDLK_RALT:
            case SDLK_LALT:
              trackball.ButtonUp(Trackball::KEY_ALT);
              break;
          }
          break;
        case SDL_MOUSEBUTTONDOWN:       
          x = event.button.x;
          y = height - event.button.y;          
          trackball.MouseDown(x, y, Trackball::BUTTON_LEFT);
        break;
        case SDL_MOUSEBUTTONUP:          
          x = event.button.x;
          y = height - event.button.y;          
          trackball.MouseUp(x, y, Trackball::BUTTON_LEFT);
          break;
        case SDL_MOUSEMOTION: 
          while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
          x = event.motion.x;
          y = height - event.motion.y;
          trackball.MouseMove(x, y);        
          break;  
      }
        
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 0.1, 100);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0,0,6,   0,0,0,   0,1,0);        
    glRotatef(130, 1, 1, 0);
    glTranslatef(0, 1, 1);
    
    
    //    trackball.SetPosition(Similarityf(Point3f(1, 0, 0)));
    //    trackball.local.sca = 0.5;
    trackball.GetView();
    trackball.Apply();
    trackball.Draw();

    glColor3f(0, 1, 0);
    glutWireCube(1);
    
    SDL_GL_SwapBuffers();
  }

        
  // Clean up
  SDL_Quit();
  

  return -1;
}
  
