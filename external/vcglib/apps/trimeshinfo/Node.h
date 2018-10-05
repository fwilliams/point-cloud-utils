
#ifndef NODE_H
#define NODE_H
#include <list>


using namespace std;

class Node
{
public:
	//int node_type;
	virtual ~Node(void){};
	virtual void printNode()=0;
	virtual int qualifyNode()=0;
};


class NodeGroup :public Node
{ 
public:
	virtual ~NodeGroup();
	typedef list<Node *>::iterator iterator;
	list<Node *> Sons;
	virtual void addNode(Node* nd);
	virtual void printNode();
	virtual int qualifyNode();
	
};

NodeGroup::~NodeGroup() //distruttore: disalloca tutti i figli
{
	//for(iterator i=Sons.begin();i!=Sons.end();++i)
	//	delete (*i);
}

void NodeGroup::addNode(Node* nd)
{
	Sons.push_back(nd);
}
void NodeGroup::printNode()
{}

int NodeGroup::qualifyNode()
{return 0;}

const int		MAIN_NODE=					0;
const int		SLOTS_NODE=				1; 
const int		SLOT_NODE=					2; 
const int		OWNSLOT_NODE=			3; 
const int		ENTRY_NODE=				4; 
const int		VALUE_NODE=				5; 

const int		CLASSES_NODE=			6; 
const int		CLASS_NODE=				7; 
const int		OWNSLOTS_NODE=			8; 

const int		INSTANCES_NODE=		9;
const int		INSTANCE_NODE=			10;


enum values {VALUE_INTEGER, VALUE_FLOAT, VALUE_BOOL, VALUE_STRING};

//#define		MAIN_NODE					0; 
//#define		SLOTS_NODE				1; 
//#define		SLOT_NODE					2; 
//#define		OWNSLOT_NODE			3; 
//#define		ENTRY_NODE				4; 
//#define		VALUE_NODE				5; 
//
//#define		CLASSES_NODE			6; 
//#define		CLASS_NODE				7; 
//#define		OWNSLOTS_NODE			8; 
//
//#define		INSTANCES_NODE		9;
//#define		INSTANCE_NODE			10;

#endif
