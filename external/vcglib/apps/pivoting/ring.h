#ifndef RING_H
#define RING_H
// Making a "ring" data structure from the STL
#include <list>
using namespace std;

template<class T>
class ring {
  list<T> lst;
public:

  ring &operator=(const ring *r) {
     lst = r->lst;
  }
  // Declaration necessary so the following 
  // 'friend' statement sees this 'iterator' 
  // instead of std::iterator:
  class iterator;
  friend class iterator;
  class iterator : 
     public std::iterator<std::bidirectional_iterator_tag, T, ptrdiff_t> {

   public:
    typename list<T>::iterator it;
    list<T>* r;

    
    // "typename" necessary to resolve nesting:

    iterator(): r(NULL) {}
    iterator(list<T>& lst, const typename list<T>::iterator& i)
        : r(&lst), it(i) {}

    iterator &operator=(const iterator& x) {
       it = x.it;
       r = x.r;
       return *this;
    }

    bool operator==(const iterator& x) const {
      return it == x.it;
    }

    bool operator!=(const iterator& x) const {
      return !(*this == x);
    }

    typename list<T>::reference operator*() const {
      return *it;
    }

    iterator& operator++() {
      ++it;
      if(it == r->end())
        it = r->begin();
      return *this;
    }

    iterator operator++(int) {
      iterator tmp = *this;
      ++*this;
      return tmp;
    }

    iterator& operator--() {
      if(it == r->begin())
        it = r->end();
      --it;
      return *this;
    }

    iterator operator--(int) {
      iterator tmp = *this;
      --*this; 
      return tmp;
    }

    iterator operator+(int i) {
      iterator tmp = *this;
      while(i--) ++tmp;
      return tmp;
    }

    iterator operator-(int i) {
      iterator tmp = *this;
      while(i--) --tmp;
      return tmp;
    }

    iterator insert(const T& x){
      return iterator(*r, r->insert(it, x));
    }
    
    iterator erase() {
      iterator tmp = iterator(*r, r->erase(it));
      
      if (tmp.it == r->end()) tmp.it = r->begin();

      return tmp;
    }
  };
  
  iterator insert(iterator &i, const T& x) {
    typename list<T>::iterator it;
    it = lst.insert(i.it, x);
    return iterator(lst, it);
  }

  iterator push_back(const T& x) {
    typename list<T>::iterator it;
    it = lst.insert(lst.end(), x);
    return iterator(lst, it);
    
  }
 /* void push_back(const T& x) {
    lst.push_back(x);
  }*/
  
  iterator begin() {
    return iterator(lst, lst.begin());
  }
  
  int size() { return lst.size(); }
  void clear() { lst.clear(); }
  void reverse() { lst.reverse(); }
  void erase(iterator &i) { lst.erase(i.it); }
};

#endif
