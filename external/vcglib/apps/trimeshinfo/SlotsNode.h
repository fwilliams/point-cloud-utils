
#include <utility>
#include <map>
#include <list>
#include "Node.h"

using namespace std;

class ValueNode: public Node
{
public:
	ValueNode(void){node_type = VALUE_NODE;value = "empty";};
	int node_type;
	string value; //tra due tag
	virtual void printNode();
	virtual int qualifyNode();

	void setValue(ValueNode vn){value = vn.value;};
	void setValue(const char* cvn){value = cvn;};
};

void ValueNode::printNode()
{
	cout<<"ValueNode: Node type is "<<node_type<<"\n";
	cout<<"ValueNode: Node value is "<<value<<"\n";
}

int ValueNode::qualifyNode()
{return node_type;}

class EntryNode: public Node
{
public:
	EntryNode(void){node_type = ENTRY_NODE; type = "empty";};
	int node_type;
	const char* type;
	ValueNode value;
	void addValue(ValueNode vn);
	void setEntry(EntryNode en);
	virtual void printNode();
	virtual int qualifyNode();
};
void EntryNode::addValue(ValueNode vn)
{
	value.setValue(vn);
}
void EntryNode::setEntry(EntryNode en)
{
	type = en.type;
	addValue(en.value);
}

void EntryNode::printNode()
{
	cout<<"EntryNode: Node type is "<<node_type<<"\n";
	cout<<"EntryNode: Node attr. type is "<<type<<"\n";
	value.printNode();
}

int EntryNode::qualifyNode()
{return node_type;}

class OwnSlotNode: public Node
{
public:
	OwnSlotNode(void){node_type = OWNSLOT_NODE; name = "empty";};
	int node_type;
	const char* name;
	EntryNode entry;
	virtual void printNode();
	virtual int qualifyNode();
	void addEntry(EntryNode en);
	void setName(const char* s){name = s;};
};

void OwnSlotNode::printNode()
{
	cout<<"OwnSlotNode: Node type is "<<node_type<<"\n";
	cout<<"OwnSlotNode: Node name is "<<name<<"\n";
	entry.printNode();
}

int OwnSlotNode::qualifyNode()
{return node_type;}

void OwnSlotNode::addEntry(EntryNode en)
{
	entry.setEntry(en);
}

class SlotNode: public Node
{
public:
	SlotNode(void){node_type = SLOT_NODE;};
	int node_type;
	NodeGroup	own_slot;
	virtual void printNode();
	virtual int qualifyNode();
	void addOwnSlot(OwnSlotNode* os);

};

void SlotNode::addOwnSlot(OwnSlotNode* os)
{
	//OwnSlotNode* osn = new OwnSlotNode;
//	own_slots.Sons.push_back(new OwnSlotNode);
	
		own_slot.addNode(os);
}

void SlotNode::printNode()
{
	cout<<"SlotNode: Node type is "<<node_type<<"\n";
	list<Node*>::iterator it;
	for(it = own_slot.Sons.begin(); it!=own_slot.Sons.end(); ++it)
		(*it)->printNode();
}

int SlotNode::qualifyNode()
{return node_type;}

class SlotsNode: public Node
{
public:
	SlotsNode(void){node_type = SLOTS_NODE;};
	int node_type;
	NodeGroup slot;
	void addSlot(SlotNode* sn);
	virtual void printNode();
	virtual int qualifyNode();
};

void SlotsNode::addSlot(SlotNode* sn)
{
	slot.Sons.push_back(new SlotNode);
	SlotNode* slp = (SlotNode*) slot.Sons.front();
	list<Node*>::iterator it;
	for(it = sn->own_slot.Sons.begin(); it!=sn->own_slot.Sons.end(); ++it)
		slp->addOwnSlot(((OwnSlotNode*)(*it)));
}

void SlotsNode::printNode()
{
	cout<<"SlotsNode: Node type is "<<node_type<<"\n";
	list<Node*>::iterator it;
	for(it = slot.Sons.begin(); it!=slot.Sons.end(); ++it)
		(*it)->printNode();
}

int SlotsNode::qualifyNode()
{return node_type;}
