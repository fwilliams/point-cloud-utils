/* 

This extends the example in sample/trimesh_sdl.h.
It loads the mesh specified in the command line and shows it. 

Navigate around the object with WASD (doom-like) or arrow modes.
(wheel, or PGUP, PGDOWN, moves vertically)

[F1] changes rendering style.
[SPACE] toggles trackball modes.


This example shows:
  
- the use of a trackball in
a time-syncronized application like a game.
(exactly FPS times every seconds, the trackball is animated,
but rendering can occurs fewer times).

- how to set a trackball
with NavigatorAwsd mode, i.e. to NAVIGATE around or inside
the object rather than rotate it in front of the camera.

*/

#include <SDL/SDL.h>
#include <gl/glew.h>
#include <vector>

/*include the base definition for the vertex, face, and meshes */
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/complex/complex.h>

/*include the algorihm that update bounding box and normals*/
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>

/*include the importer from disk*/
#include <wrap/io_trimesh/import.h>

/*include wrapping of the trimesh towards opengl*/
#include <wrap/gl/trimesh.h>

/*include the trackball: NOTE, you the implementation of the trackball is found in the files:
wrap/gui/trackball.cpp and wrap/gui/trackmode.cpp. You should include these files in your solution
otherwise you'll get linking errors */
#include <wrap/gui/trackball.h>

/* we define our own SDL_TIMER event */
#define SDL_TIMER SDL_USEREVENT

// How many FPS?
const int FPS = 30; 

using namespace vcg;
using namespace std;

/* Definition of our own types for Vertices, Faces, Meshes*/
class CEdge;    // dummy prototype never used
class CFace;
class CVertex : public VertexSimp2< CVertex, CEdge, CFace, vertex::Coord3f, vertex::Normal3f >{};
class CFace   : public FaceSimp2<   CVertex, CEdge, CFace, face::VertexRef, face::Normal3f > {};
class CMesh   : public vcg::tri::TriMesh< vector<CVertex>, vector<CFace> > {};



////////////////////////////////////////////////////////////////////////////
// Globals: the mesh, the OpenGL wrapper to draw the mesh and the trackball.
CMesh mesh;
vcg::GlTrimesh<CMesh> glWrap;
vcg::Trackball track;

/* Internal state of the application */
int drawMode=2;
int trackballMode=1;
int width =800;
int height = 600;
vcg::Point3f observerPos(0,0.2,3); // (initial) point of view (object is in 0,0,0)
vcg::Point3f objectPos(0,0,0); // (initial) point of view (object is in 0,0,0)


/* Simple function that renders a floor */
void RenderFloor(){
  glHint(GL_FOG_HINT, GL_NICEST);	
  glDisable(GL_LIGHTING);
  double K = 500, c=0.7;
  glBegin(GL_QUADS);
  glColor3f(c,c,c);
  glNormal3f(0,1,0);
  glVertex3d(+K,0,+K);
  glVertex3d(+K,0,-K);
  glVertex3d(-K,0,-K);
  glVertex3d(-K,0,+K);
  glEnd();
  glHint(GL_FOG_HINT, GL_FASTEST);	
  glEnable(GL_LIGHTING);
}

/* Sets the trackball in a new mode: 
  standard mode (rotate object in front of camera) 
  or Navigation mode (navigate around/inside object)
*/
void SetTrackball(int mode){
  // we define all possible trackModes that we could be using (static)
	static vcg::PolarMode polarMode;
	static vcg::SphereMode sphereMode;
	static vcg::NavigatorWasdMode navigatorMode;
	static vcg::InactiveMode inactiveMode;
	static vcg::ScaleMode scaleMode;
	static vcg::PanMode panMode;	
	static vcg::ZMode zMode;	
	
	// customize navigation mode... 
	navigatorMode.SetTopSpeedsAndAcc(1.2f,0.6f,6);
	// this adds a neat human stepping effect
	navigatorMode.SetStepOnWalk(0.5f,0.015f);
	
	track.modes.clear();
	track.Reset();
	
  switch (mode) {
  case  1:
    // switch to navigator mode
	  track.modes[vcg::Trackball::BUTTON_NONE] = NULL;
	  track.modes[vcg::Trackball::WHEEL] = 
	  track.modes[vcg::Trackball::BUTTON_LEFT] = 
	  track.idle_and_keys_mode = &navigatorMode;

    track.inactive_mode = NULL; // nothing to draw when inactive
    track.SetPosition( observerPos );
    break;
  default:
    // sweitch to default trackmode -- this is equivalent to a call to track->SetDefault()
    track.modes[vcg::Trackball::BUTTON_NONE] = NULL;
    track.modes[vcg::Trackball::BUTTON_LEFT] = &sphereMode;
    track.modes[vcg::Trackball::BUTTON_LEFT | vcg::Trackball::KEY_CTRL] = 
    track.modes[vcg::Trackball::BUTTON_MIDDLE] = &panMode;
    track.modes[vcg::Trackball::WHEEL] = 
    track.modes[vcg::Trackball::BUTTON_LEFT | vcg::Trackball::KEY_SHIFT] = &scaleMode;
    track.modes[vcg::Trackball::BUTTON_LEFT | vcg::Trackball::KEY_ALT] = &zMode;
    track.modes[vcg::Trackball::BUTTON_MIDDLE | vcg::Trackball::KEY_ALT] = 
	  track.idle_and_keys_mode = &sphereMode;
	  
    track.inactive_mode = &inactiveMode; // draw a sphere when inactive
    track.SetPosition( objectPos );
  }
}

/* Sets up OpenGL state */
void initGL()
{
  glClearColor(0, 0, 0, 0); 
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_NORMALIZE);
  glEnable(GL_FOG);
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT);
  glDisable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
}

/* Response to a window resize event. */
void myReshapeFunc(GLsizei w, GLsizei h)
{
  SDL_SetVideoMode(w, h, 0, SDL_OPENGL|SDL_RESIZABLE|SDL_DOUBLEBUF);
	width=w;
  height=h;
  glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
  initGL();
}

/* sends a SDL_TIMER event to self (timer callback) */
Uint32 timerCallback(Uint32 interval, void *param){ 
  static SDL_Event e; // the msg
  e.type = SDL_TIMER; // its content
  SDL_PushEvent(&e);
  return interval;
}
 
/* sends a redraw event to self */
void sendRedraw() {
  static SDL_Event e; // the msg
  e.type = SDL_VIDEOEXPOSE; // its content
  SDL_PushEvent(&e);
}

/* clears any "redraw" even still in the event queue */
void drainRedrawEvents(){
  static SDL_Event tmp[200];
  int eventRemoved = SDL_PeepEvents(tmp,200,SDL_GETEVENT, SDL_EVENTMASK(SDL_VIDEOEXPOSE));
}

/* Response to a redraw event: renders the scene */
void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, width/(float)height, 0.01, 10);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslate( -observerPos );

		track.GetView();
	track.Apply();
    glPushMatrix();
    float d=mesh.bbox.Diag();
    glScale(1.5f/d);
    Point3f p = glWrap.m->bbox.Center();
    p[1] = glWrap.m->bbox.min[1];
		glTranslate(-p);	

		// the trimesh drawing calls
		switch(drawMode)
		{
		  case 0: glWrap.Draw<vcg::GLW::DMSmooth,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 1: glWrap.Draw<vcg::GLW::DMWire,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 2: glWrap.Draw<vcg::GLW::DMFlatWire, vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 3: glWrap.Draw<vcg::GLW::DMHidden,   vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  case 4: glWrap.Draw<vcg::GLW::DMFlat,     vcg::GLW::CMNone,vcg::GLW::TMNone> ();break;
		  default: break;
		}
    glPopMatrix();
    RenderFloor();
    track.DrawPostApply();
    SDL_GL_SwapBuffers();
}

/* Response to a timer event: animates the status of the game (and the trackball) */
// - retunrs true if anything changed.
bool onTimer(){
  int res = false;
  if ( track.IsAnimating() ) {
    track.Animate(1000/FPS);
    res = true;
  }
  // insert any other animation processing here
  return res;
}


/* Helper function: translates SDL codes into VCG codes */
vcg::Trackball::Button SDL2VCG(int code){
  switch (code) {
    case SDL_BUTTON_LEFT: return vcg::Trackball::BUTTON_LEFT;
    case SDL_BUTTON_MIDDLE: return vcg::Trackball::BUTTON_MIDDLE;
    case SDL_BUTTON_RIGHT: return vcg::Trackball::BUTTON_RIGHT;
		case SDLK_RCTRL:
		case SDLK_LCTRL:  return vcg::Trackball::KEY_CTRL;
		case SDLK_RALT:
		case SDLK_LALT:  return vcg::Trackball::KEY_ALT;
		case SDLK_LSHIFT:
		case SDLK_RSHIFT:  return vcg::Trackball::KEY_SHIFT; 
		case SDLK_LEFT:
    case SDLK_a:  return vcg::Trackball::KEY_LEFT;
		case SDLK_RIGHT:
    case SDLK_d:  return vcg::Trackball::KEY_RIGHT;
		case SDLK_UP:
    case SDLK_w:  return vcg::Trackball::KEY_UP;
		case SDLK_DOWN:
    case SDLK_s:  return vcg::Trackball::KEY_DOWN;
		case SDLK_PAGEUP:
    case SDLK_r:  return vcg::Trackball::KEY_PGUP;
		case SDLK_PAGEDOWN:
    case SDLK_f: return vcg::Trackball::KEY_PGDOWN;
    return vcg::Trackball::BUTTON_NONE;
  }
}

/* initializes SDL */
bool initSDL(const std::string &str) {
  if ( SDL_Init(SDL_INIT_VIDEO|SDL_INIT_TIMER) < 0 ) {
    fprintf(stderr, "Unable to initialize SDL: %s\n", SDL_GetError());
    exit(1);
  }
  if ( SDL_SetVideoMode(width, height, 0, SDL_OPENGL|SDL_RESIZABLE) == NULL ) {
    fprintf(stderr, "Unable to create OpenGL screen: %s\n", SDL_GetError());
    SDL_Quit();
    exit(2);
  }
  
  SDL_AddTimer(1000/FPS, timerCallback, NULL); 
  
  SDL_WM_SetCaption(str.c_str(), str.c_str());  
  myReshapeFunc(width, height);
  return true;
}

/* The main event Loop */
int sdlLoop() {
  bool quit=false;

  bool redraw_needed = false; // true whan a scene needs a redraw
  
 	SDL_Event event;
	while( !quit ) {

    SDL_WaitEvent(&event);
    switch( event.type ) {
      case SDL_QUIT:  quit = true; break; 
      case SDL_VIDEORESIZE : myReshapeFunc(event.resize.w,event.resize.h); break;
      case SDL_KEYDOWN: 
        switch (event.key.keysym.sym) {
			    case SDLK_ESCAPE: exit(0); break;	
			    case SDLK_F1: drawMode= (drawMode+1)%5; printf("Current Mode %i\n",drawMode); break;	
			    case SDLK_HOME: track.Reset(); break;
          case SDLK_SPACE: 
            trackballMode= (trackballMode+1)%2; printf("Trackball Mode %i\n",drawMode); 
            SetTrackball(trackballMode);
            break;	
          default: track.ButtonDown( SDL2VCG( event.key.keysym.sym) );
			  }
			  redraw_needed = true;
        break;
      case SDL_KEYUP: 
			  track.ButtonUp( SDL2VCG( event.key.keysym.sym) ); break;
      case SDL_MOUSEBUTTONDOWN:   
	      switch(event.button.button) {
          case SDL_BUTTON_WHEELUP:    track.MouseWheel( 1); redraw_needed = true; break;
          case SDL_BUTTON_WHEELDOWN:  track.MouseWheel(-1); redraw_needed = true; break;
          default: track.MouseDown(event.button.x, (height - event.button.y), SDL2VCG(event.button.button) ); break;
        } break;
      case SDL_MOUSEBUTTONUP:   
        track.MouseUp(event.button.x, (height - event.button.y), SDL2VCG(event.button.button) ); break;       
        break;
      case SDL_MOUSEMOTION: 
	      while(SDL_PeepEvents(&event, 1, SDL_GETEVENT, SDL_MOUSEMOTIONMASK));
	      track.MouseMove(event.button.x, (height - event.button.y));
 			  redraw_needed = true;
	      break;
	    case SDL_TIMER:
        if (onTimer())  redraw_needed = true;
        if (redraw_needed) sendRedraw(); // justs sends a redraw event!
 			  redraw_needed = false;
        break;
      case SDL_VIDEOEXPOSE:
        // any rendering is done ONLY here.
		    display();
		    drainRedrawEvents();
		    break;
      default: break;
      }

	}

  SDL_Quit();
  return -1;
}


int main(int argc, char *argv[]) {	
	// Generic loading of the mesh from disk
	if(vcg::tri::io::Importer<CMesh>::Open(mesh,argv[1])!=0) {
      fprintf(stderr,"Error reading file %s\n",argv[1]);
			exit(0);
		}

  // Initialize the mesh itself
	vcg::tri::UpdateBounding<CMesh>::Box(mesh);      // update bounding box
  //vcg::tri::UpdateNormals<CMesh>::PerVertexNormalizePerFaceNormalized(mesh); // update Normals
  vcg::tri::UpdateNormals<CMesh>::PerVertexPerFace(mesh); // update Normals
	
	// Initialize the wrapper
  glWrap.m = &mesh;
  glWrap.Update();
  
  SetTrackball(trackballMode);
  
  // we will do exaclty an animation step every 1000/FPS msecs.
  track.SetFixedTimesteps(true);
	
  initSDL("SDL_minimal_viewer");
  initGL();
  sdlLoop();
	exit(0);
}



