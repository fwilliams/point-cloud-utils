#ifndef VCG_PIVOT_H
#define VCG_PIVOT_H

#include <vector>
#include <list>


#include "vcg/space/index/grid_static_ptr.h"
#include "vcg/complex/algorithms/closest.h"

namespace vcg {
namespace tri {


struct Hinge { 
      int v0, v1, v2;   //v0, v1 represent the Hinge, v2 the other vertex in the face
                        //this Hinge belongs to
      int face;        //corresponding face
      Point3f center;  //center of the sphere touching the face
      int count;   //test delay touch Hinges.  
      
      //the loops in the front are mantained as a double linked list
      std::list<Hinge>::iterator next;            
      std::list<Hinge>::iterator previous;
    
      Hinge() {}
      Hinge(int _v0, int _v1, int _v2, int _face, Point3f &_center): 
               v0(_v0), v1(_v1), v2(_v2), 
               face(_face), center(_center), count(0) {
                 assert(v0 != v1 && v1 != v2 && v0 != v2);
               }
    };
    
template <class MESH>
class Pivot {
  public:
//    typedef CMesh MESH;
    typedef GridStaticPtr<typename MESH::VertexType, typename MESH::ScalarType > StaticGrid;
    
    

    float radius;  //default 1 (not meaningful
    float mindist; //minimum distance between points in the mesh (% of radius)   
    float crease;  // -0.5 
    Box3f box;
    
    MESH &mesh;    
    StaticGrid grid;
    
    /* front Hinges of the mesh: 
       to expand the front we get the first Hinge
       if an Hinge cannot create a new triangle it is marked dead and moved
       to the end of the list
       the new Hinges are inserted to the back (before dead_begin) */
       
    std::list<Hinge> front;   
    std::list<Hinge> deads;
    std::vector<int> nb; //number of fronts a vertex is into,
                         //this is used for the Visited and Border flags
                         //but adding topology may not be needed anymode

      
    Pivot(MESH &_mesh, float _radius, float _mindist = 0.05, float _crease = -0.5): 
       mesh(_mesh), radius(_radius), mindist(_mindist), crease(_crease) {
    
      //Compute bounding box. (this may be passed as a parameter?
      for(int i = 0; i < mesh.vert.size(); i++)
        box.Add(mesh.vert[i].P());
     
      /* we need to enlarge the grid to allow queries from little outside of the box 
         Someone is a bit lazy... */        
      box.Offset(4*radius);
      grid.Set(mesh.vert.begin(), mesh.vert.end(), box);
      nb.clear();
      nb.resize(mesh.vert.size(), 0);
      for(int i = 0; i < mesh.vert.size(); i++) 
         mesh.vert[i].ClearFlags();
      
    }
    
    
    /* select a vertex at random, a small group of nearby vertices and looks
       for a sphere that touches 3 and contains none.
       Use the center of the box to get a sphere inside (or outside) the model 
       You may be unlucky... */
       
    bool seed(bool outside = true, int start = -1) {         
      //pick a random point (well...)
      if(start == -1) start = rand()%mesh.vert.size();
      
      //get a sphere of neighbours
      std::vector<int> targets;
      std::vector<float> dists;
      int n = getInSphere(mesh.vert[start].P(), 2*radius, targets, dists);
      if(n < 3) { 
        //bad luck. we should call seed again (assuming random pick) up to
        //some maximum tries. im lazy.
        return false;
      }
      int v0, v1, v2;
      bool found = false;
      //find a triplet that does not contains any other point
      Point3f center;
      for(int i = 0; i < n; i++) {
        v0 = targets[i];
        CVertex &vv0 = mesh.vert[v0];
        if(vv0.IsD() || vv0.IsB() || vv0.IsV()) continue;
        Point3f &p0 = mesh.vert[v0].P();
        Point3f out = (p0 - box.Center());
        if(!outside) out = -out;
    
        for(int k = i+1; k < n; k++) {
          v1 = targets[k];            
          CVertex &vv1 = mesh.vert[v1];
          if(vv1.IsD() || vv1.IsB() || vv1.IsV()) continue;
          Point3f &p1 = mesh.vert[v1].P();            
          if((p1 - p0).Norm() < mindist*radius) continue;
    
          for(int j = k+1; j < n; j++) {
            v2 = targets[j];
            CVertex &vv2 = mesh.vert[v2];
            if(vv2.IsD() || vv2.IsB() || vv2.IsV()) continue;
            Point3f &p2 = mesh.vert[v2].P();
            if((p2 - p0).Norm() < mindist*radius) continue;
            if((p2 - p1).Norm() < mindist*radius) continue;
            Point3f normal = (p1 - p0)^(p2 - p0);
            //check normal pointing inside
            if(normal * out < 0) continue;
            if(!findSphere(p0, p1, p2, center)) continue;
            
            bool failed = false;
            //check no other point inside
            for(int t = 0; t < n; t++) {
              Point3f &p = mesh.vert[targets[t]].P();
              if((center - p).Norm() <= radius) {
                failed = true;
                break;
              }
            }
            if(failed) continue;  
            found = true;
            i = k = j = n;
          }
        }
      }
      
      if(!found)  //see bad luck above
        return false;
      
      assert(!front.size());
      //TODO: should i check the Hinges too?
      addFace(v0, v1, v2);
      
      //create the border of the first face  
      std::list<Hinge>::iterator e = front.end();
      std::list<Hinge>::iterator last;
      for(int i = 0; i < 3; i++) {
        int v0 = (int)(mesh.face.back().V0(i));
        int v1 = (int)(mesh.face.back().V1(i));    
        int v2 = (int)(mesh.face.back().V2(i));    
        nb[v0] = 1;
        assert(!mesh.vert[v0].IsB());
        mesh.vert[v0].SetB();
        Hinge Hinge(v0, v1, v2, 0, center);
        Hinge.previous = e;
        e = front.insert(front.begin(), Hinge);
        if(i == 0) last = e;
        (*Hinge.previous).next = e;
        
        cluster(v0);
      } 
    
      //connect last and first
      (*e).next = last;
      (*last).previous = e;
      return true;
    }
      
         
    /* expand the front adding 1 face. Return false on failure (id when
       all Hinges are dead  returns:
       1: added a face
       0: added nothing
       -1: finished */
    int addFace() {
      //We try to seed again
      if(!mesh.face.size()) {
        for(int i = 0; i < 100; i++) 
          if(seed()) return 1;
        return -1;
      }
      
      if(!front.size()) {
        //maybe there are unconnected parts of the mesh:
        //find a non D, V, B point and try to seed if failed D it.
        for(int i = 0; i < mesh.vert.size();i ++) {
          CVertex &v = mesh.vert[i];
          if(v.IsD() || v.IsV() || v.IsB()) continue;
          if(seed(true, i)) return 1;
          v.SetD();
        }
        return -1;
      }
    
      std::list<Hinge>::iterator ei = front.begin();
      Hinge &e = *ei;
      Hinge &previous = *e.previous;           
      Hinge &next = *e.next;  
      int v0 = e.v0, v1 = e.v1;
    
      //last triangle missing. or it is the first?
      if(0 &&next.next == e.previous) {  
    
        int v[3] = { previous.v0, next.v0, e.v0 };
        int c[3] = { 0, 0, 0 };
            
        for(int k = 0; k < 3; k++) {
          int vert = v[k];
          nb[vert]--;
          if(nb[vert] == 0) {       
            mesh.vert[vert].SetV();
            mesh.vert[vert].ClearB();
          }              
        }        
        assert(previous.previous == e.next); 
        addFace(previous.v0, next.v0, e.v0);                       
    
        front.erase(e.previous);
        front.erase(e.next);        
        front.erase(ei);    
        
        return 1;
      }
    
      int v2;
      Point3f center;
    
      std::vector<int> targets;
      bool success = pivot(e, v2, center, targets);
    
      //if no pivoting move this thing to the end and try again
      //or we are trying to connect to the inside of the mesh. BAD.
      if(!success || mesh.vert[v2].IsV()) {
        killHinge(ei);
        return 0;
      } 
    
      //does v2 belongs to a front? (and which?)
      std::list<Hinge>::iterator touch = touches(v2, ei);
    
      assert(v2 != v0 && v2 != v1);  
    
      int fn = mesh.face.size();
      if(touch != front.end()) {       
    
        if(!checkHinge(v0, v2) || !checkHinge(v2, v1)) {                      
          killHinge(ei);
          return 0;
        }
        
        if(v2 == previous.v0) {    
          /*touching previous Hinge  (we reuse previous)        
                                    next
             ------->v2 -----> v1------>
                      \       /
                       \     /
               previous \   / e
                         \ /
                          v0           */
          
          detach(v0);
        
          previous.v1 = v1;       
          previous.v2 = v0;
          previous.face = fn;
          previous.center = center;
          
          previous.next = e.next;
          next.previous = e.previous;
          moveBack(e.previous);            
          
          //this checks if we can glue something to e.previous
          trovamiunnome(e.previous);      
          front.erase(ei);    
          
           
        } else if(v2 == next.v1) {
        /*touching previous Hinge  (we reuse next)        
          previous
             ------->v0 -----> v2------>
                      \       /
                       \     /
                        \   / next
                         \ /
                          v1           */      
    
          detach(v1);
          
          next.v0 = v0;
          next.v2 = v1;
          next.face = fn;
          next.center = center;
          next.previous = e.previous;
          previous.next = e.next;
    //      moveBack(e.next);
          
          //this checks if we can glue something to e.previous
          trovamiunnome(e.next);         
          front.erase(ei);
              
        } else {
    
    /*   this code would delay the joining Hinge to avoid bad situations not used but..
         if(e.count < 2) {
            e.count++;
            moveBack(ei);
            return true;         
          }*/
    
        //touching some loop: split (or merge it is local does not matter.
        //like this 
        /*                 
                    left        right
                  <--------v2-<------
                          /|\
                         /   \
                     up /     \ down
                       /       \
                      /         V
                 ----v0 - - - > v1---------
                          e                         */           
          std::list<Hinge>::iterator left = touch;
          std::list<Hinge>::iterator right = (*touch).previous;      
          std::list<Hinge>::iterator up = ei;
          
          if(v1 == (*right).v0 || v0 == (*left).v1) {
//            cout << "Bad join.\n";
            killHinge(ei);
            return 0;
          }
          
          nb[v2]++;    
                      
          std::list<Hinge>::iterator down = newHinge(Hinge(v2, v1, v0, fn, center));      
    
          (*right).next = down;
          (*down).previous = right;
    
          (*down).next = e.next;
          next.previous = down;      
    
          (*left).previous = up;
          (*up).next = left;
          
          (*up).v2 = v1;      
          (*up).v1 = v2;
          (*up).face = fn;
          (*up).center = center;
          moveBack(ei);
        }                         
    
              
      
      } else {
        assert(!mesh.vert[v2].IsB()); //fatal error! a new point is already a border?
    
        /*  adding a new vertex
                 
                           v2
                          /|\
                         /   \
                     up /     \ down
                       /       \
                      /         V
                 ----v0 - - - > v1--------- */
        cluster(v2);        
        nb[v2]++;                 
        mesh.vert[v2].SetB();
        std::list<Hinge>::iterator down = newHinge(Hinge(v2, v1, v0, fn, center));
        (*down).previous = ei;
        (*down).next = e.next;
        next.previous = down;
        
        e.v2 = v1;    
        e.v1 = v2;
        e.face = fn;
        e.center = center;
        e.next = down; 
        moveBack(ei);
      }
      addFace(v0, v2, v1);
      return 1;
    }        
    

    
    /* return new vertex and the center of the new sphere pivoting from Hinge
       if the vertex belongs to another Hinge, touch points to it. */
    bool pivot(Hinge &Hinge, int &candidate, Point3f &end_pivot, std::vector<int> &targets) {
        Point3f v0 = mesh.vert[Hinge.v0].P();
        Point3f v1 = mesh.vert[Hinge.v1].P();  
        Point3f v2 = mesh.vert[Hinge.v2].P();  
        /* TODO why using the face normals everything goes wrong? should be
           exactly the same................................................
           Check if the e.face is correct.
           Point3f &normal = mesh.face[Hinge.face].N();
        */
    
        Point3f normal = ((v1 - v0)^(v2 - v0)).Normalize();
        
        Point3f middle = (v0 + v1)/2;    
        Point3f start_pivot = Hinge.center - middle;          
        Point3f axis = (v1 - v0);
        
        float axis_len = axis.SquaredNorm();
        if(axis_len > 4*radius*radius) return false;
        axis.Normalize();
        
        // r is the radius of the thorus of all possible spheres passing throug v0 and v1
        float r = sqrt(radius*radius - axis_len/4);
        
    
        std::vector<float> dists;    
        getInSphere(middle, r + radius, targets, dists);
    
        if(targets.size() == 0) return false; //this really would be strange but one never knows.
    
        candidate = -1;
        float minangle = 0;
        Point3f center;  //to be computed for each sample
        for(int i = 0; i < targets.size(); i++) {      
          int id = targets[i];
          
          if(id == Hinge.v0 || id == Hinge.v1 || id == Hinge.v2) continue;
    
          if(mesh.vert[id].IsD()) {
            continue;      
          }
    
    
          Point3f p = mesh.vert[id].P();
    
          /* Prevent 360 Hinges, also often reject ~ 50% points */      
          Point3f n = ((p - v0)^(v1 - v0)).Normalize();
          if(n * normal < -0.5) {
            continue;               
          }
          
    
          /* Find the sphere through v0, p, v1 (store center on end_pivot */
          if(!findSphere(v0, p, v1, center)) {
            continue;      
          }
          
          /* Angle between old center and new center */
          float alpha = angle(start_pivot, center - middle, axis);
    
          /* adding a small bias to already chosen vertices.
             doesn't solve numerical problems, but helps. */
          if(mesh.vert[id].IsB()) alpha -= 0.001;
          
          /* Sometimes alpha might be little less then M_PI while it should be 0,
             by numerical errors: happens for example pivoting 
             on the diagonal of a square. */
          
          if(alpha > 2*M_PI - 0.8) {               
            // Angle between old center and new *point* 
            //TODO is this really overshooting? shouldbe enough to alpha -= 2*M_PI
            Point3f proj = p - axis * (axis * p - axis * middle);
            float beta = angle(start_pivot, proj - middle, axis);
          
            if(alpha > beta) alpha -= 2*M_PI; 
          }
          if(candidate == -1 || alpha < minangle) {        
            candidate = id; 
            minangle = alpha;
            end_pivot = center;
          }
        }
        //found no point suitable.
        if(candidate == -1) {
          return false;
        }
           
        assert(candidate != Hinge.v0 && candidate != Hinge.v1);
        return true;
    }         
    
  private:
     //front management:
     //Add a new Hinge to the back of the queue
     std::list<Hinge>::iterator newHinge(Hinge e) {                  
       return front.insert(front.end(), e);
     }     
     //move an Hinge among the dead ones
     void killHinge(std::list<Hinge>::iterator e) {
       deads.splice(deads.end(), front, e);
     }

     //move an Hinge to the back of the queue
     void moveBack(std::list<Hinge>::iterator e) {
       front.splice(front.end(), front, e);          
     }
     
     void moveFront(std::list<Hinge>::iterator e) {
       front.splice(front.begin(), front, e);
     }
    bool checkHinge(int v0, int v1) {
      int tot = 0;
      //HACK to speed up things until i can use a seach structure
      int i = mesh.face.size() - 2*(front.size());
    //  i = 0;
      if(i < 0) i = 0;
      for(; i < mesh.face.size(); i++) { 
        CFace &f = mesh.face[i];
        for(int k = 0; k < 3; k++) {
          if(v1== (int)f.V(k) && v0 == (int)f.V((k+1)%3)) ++tot;
          else if(v0 == (int)f.V(k) && v1 == (int)f.V((k+1)%3)) { //orientation non constistent
             return false;
          }
        }
        if(tot >= 2) { //non manifold
          return false;
        }
      }
      return true;
    }        
    
               
    void Pivot::cluster(int v) {
      /* clean up too close points */
        std::vector<int> targets;
        std::vector<float> dists;    
        getInSphere(mesh.vert[v].P(), mindist*radius, targets, dists);
        
        for(int i = 0; i < targets.size(); i++) {
          int id = targets[i];
          if(id == v) continue;
    
          CVertex &v = mesh.vert[id];
          if(v.IsD() || v.IsV() || v.IsB()) continue;
          v.SetD();
        }
            
    }
    
    void Pivot::trovamiunnome(std::list<Hinge>::iterator e) {
      if(glue((*e).previous, e)) return;
      glue(e, (*e).next);
    }
    
    //glue toghether a and b (where a.next = b
    bool Pivot::glue(std::list<Hinge>::iterator a, std::list<Hinge>::iterator b) {
      if((*a).v0 != (*b).v1) return false; 
      
      std::list<Hinge>::iterator previous = (*a).previous;
      std::list<Hinge>::iterator next = (*b).next;
      (*previous).next = next;
      (*next).previous = previous;
      detach((*a).v1);
      detach((*a).v0); 
      front.erase(a);
      front.erase(b);  
      return true;
    }
    
    void Pivot::detach(int v) {
      assert(nb[v] > 0);
      if(--nb[v] == 0) {
        mesh.vert[v].SetV();
        mesh.vert[v].ClearB();      
      }
    }

          
    /* compute angle from p to q, using axis for orientation */
    float angle(Point3f p, Point3f q, Point3f &axis) {
      p.Normalize();
      q.Normalize();
      Point3f vec = p^q;
      float angle = acos(p*q);
      if(vec*axis < 0) angle = -angle;
      if(angle < 0) angle += 2*M_PI;
      return angle;
    }          
    /* add a new face. compute normals. */
    void addFace(int a, int b, int c) {
      CFace face;
      face.V(0) = (CVertex *)a;
      face.V(1) = (CVertex *)b;
      face.V(2) = (CVertex *)c;
      Point3f &p0 = mesh.vert[a].P();
      Point3f &p1 = mesh.vert[b].P();
      Point3f &p2 = mesh.vert[c].P();
      face.N() = ((p1 - p0)^(p2 - p0)).Normalize();
      
      mesh.face.push_back(face);
      mesh.fn++;
    }         
              
                             
    /* intersects segment [v0, v1] with the sphere of radius radius. */
    bool intersect(int v0, int v1, Point3f &center) {
      Point3f m =  mesh.vert[v1].P() -  mesh.vert[v0].P();
      float t = m*(center - mesh.vert[v0].P());
      if(t < 0) return false;
      if(t > m*m) return false;
      return true;
    }         


    float distance(int v0, int v1, Point3f &center) {
      Point3f m =  mesh.vert[v1].P() -  mesh.vert[v0].P();
      float t = m*(center - mesh.vert[v0].P())/(m*m);
      Point3f p = mesh.vert[v0].P() + m*t;
      return (p - center).Norm();
    }             


    /* return all point in a given ball, notice as i want the index
       of the vertices not the pointers... this may change in future */
    unsigned int getInSphere(vcg::Point3f &p, float distance, 
                             std::vector<int> &results,
                             std::vector<float> &dists) {
      std::vector<CVertex *> ptr;
      std::vector<Point3f> points;
      int n = vcg::tri::GetInSphereVertex(mesh, grid, p, distance, ptr, dists, points);
      for(int i = 0; i < ptr.size(); i++) 
        results.push_back(ptr[i] - &(mesh.vert[0]));
      return n;
    }                                                


    
    /* returns the sphere touching p0, p1, p2 of radius r such that
       the normal of the face points toward the center of the sphere */
    bool findSphere(Point3f &p0, Point3f &p1, Point3f &p2, Point3f &center) {
      Point3f q1 = p1 - p0;
      Point3f q2 = p2 - p0;  
    
      Point3f up = q1^q2;
      float uplen = up.Norm();
    
      //the three points are aligned
      if(uplen < 0.001*q1.Norm()*q2.Norm()) return false;
      up /= uplen;
      
    
      float a11 = q1*q1;
      float a12 = q1*q2;
      float a22 = q2*q2;
    
      float m = 4*(a11*a22 - a12*a12);
      float l1 = 2*(a11*a22 - a22*a12)/m;
      float l2 = 2*(a11*a22 - a12*a11)/m;
    
      center = q1*l1 + q2*l2;
      float circle_r = center.Norm();
      if(circle_r > radius) return false; //need too big a sphere
    
      float height = sqrt(radius*radius - circle_r*circle_r);
      center += p0 + up*height;
    
      return true;
    }         
    std::list<Hinge>::iterator touches(int v, std::list<Hinge>::iterator e) {
      //TODO what happens when it touches more than one front?
      //might still work.
    
      std::list<Hinge>::iterator touch = front.end();
      if(mesh.vert[v].IsB()) {
        //test nearby Hinges: it is faster
        std::list<Hinge>::iterator p = e;
        p = (*e).previous;
        if(v == (*p).v0) return p;
        e = (*e).next;
        if(v == (*e).v0) return e;
    
        p = (*p).previous;
        if(v == (*p).v0) return p;
        e = (*e).next;
        if(v == (*e).v0) return e;
    
        //test all. sigh.
    
        for(std::list<Hinge>::iterator k = front.begin(); k != front.end(); k++) {
          if(v == (*k).v0) { 
            touch = k;
            break;
          }
        }
        for(std::list<Hinge>::iterator k = deads.begin(); k != deads.end(); k++) {
          if(v == (*k).v0) { 
            touch = k;
            break;
          }
        }
        assert(touch != front.end());
      }
    
      return touch;
    }                              
    
    
 public:
        
 };
        

}//namespace
}//namespace
/*  CODE FOR PIVOTING IN  A TOUCH SITUATION not used now.

    //if touch we want to check the ball could really pivot around that point
    if(touch != front.end() && touch != (*Hinge.next).next && touch != Hinge.previous) {
      Point3f &hinge = mesh.vert[min].P();      
      Point3f target = (*touch).center - hinge;
      float d = (target * start_pivot)/(target.Norm()*start_pivot.Norm());
      if(d < -0.8) {
        return false;
      }
      
      

      if(d < 0.5) { //they are far enough so test .
        Point3f naxis = (target ^ start_pivot).Normalize();
        Point3f d1 = naxis^start_pivot;
        Point3f d2 = target^naxis;          
        
        for(int i = 0; i < targets.size(); i++) {
          int id = targets[i];
          if(id == Hinge.v0 || id == Hinge.v1 || id == Hinge.v2 || id == min) continue;
          if(mesh.vert[id].IsD()) {
            continue;
          }

          Point3f intruder = mesh.vert[targets[i]].P() - hinge;
          float h = intruder*naxis;
          if(fabs(h) > radius) continue;
          intruder -= naxis*h;
          assert(fabs(intruder *naxis) < 0.01);
          float off = radius - intruder.Norm(); //(distance from the center ring of the thorus
          if(h*h + off*off > radius*radius) continue; //outside of thorus
          if(d1*intruder < 0 || d2*intruder < 0) continue; //ouside of sector
          cout << "could not pivot while touching;\n";
          return false;
        }
        
      }          
    }*/


#endif
