#include <utility>
#include <map>
#include <list>
#include "Node.h"

using namespace std;


class OwnSlotsNode: public Node
{
public:
	OwnSlotsNode(void){node_type = OWNSLOTS_NODE;};
	int node_type;
	NodeGroup	own_slot;
	virtual void printNode();
	virtual int qualifyNode();
	void addOwnSlot(NodeGroup* ng);
	void addOwnSlot(OwnSlotNode* os);
};

void OwnSlotsNode::addOwnSlot(NodeGroup* ng)
{	
	list<Node*>::iterator it;
	for(it = ng->Sons.begin(); it!=ng->Sons.end(); ++it)
		own_slot.addNode(*it);
}

void OwnSlotsNode::addOwnSlot(OwnSlotNode* os)
{
	//OwnSlotNode* osn = new OwnSlotNode;
//	own_slots.Sons.push_back(new OwnSlotNode);
	
		own_slot.addNode(os);
}


void OwnSlotsNode::printNode()
{
	cout<<"OwnSlotsNode: Node type is "<<node_type<<"\n";
	list<Node*>::iterator it;
	for(it = own_slot.Sons.begin(); it!=own_slot.Sons.end(); ++it)
		(*it)->printNode();
}

int OwnSlotsNode::qualifyNode()
{return node_type;}

class ClassNode: public Node
{
public:
	ClassNode(void){node_type = CLASS_NODE;};
	int node_type;
	OwnSlotsNode own_slots;

	void addOwnSlots(OwnSlotsNode* own_slots);
virtual void printNode();
virtual int qualifyNode();
};

void ClassNode::addOwnSlots(OwnSlotsNode* sn)
{
	own_slots.addOwnSlot(&sn->own_slot);
}

void ClassNode::printNode()
{
	cout<<"ClassNode: Node type is "<<node_type<<"\n";
	own_slots.printNode();
}

int ClassNode::qualifyNode()
{return node_type;}

class ClassesNode: public Node
{
	public:
	ClassesNode(void){node_type = CLASSES_NODE;};
	int node_type;
	NodeGroup classn;
	void addClass(ClassNode* cn);
	virtual void printNode();
	virtual int qualifyNode();
};

void ClassesNode::addClass(ClassNode* cn)
{
	classn.Sons.push_back(new ClassNode);
	ClassNode* clp = (ClassNode*) classn.Sons.front();

	clp->addOwnSlots(&cn->own_slots);


}

void ClassesNode::printNode()
{
	cout<<"ClassesNode: Node type is "<<node_type<<"\n";
	list<Node*>::iterator it;
	for(it = classn.Sons.begin(); it!=classn.Sons.end(); ++it)
		(*it)->printNode();
	
}

int ClassesNode::qualifyNode()
{return node_type;}



