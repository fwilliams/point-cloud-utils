
#include <utility>
#include <map>
#include <list>
#include "Node.h"


class InstanceNode: public Node
{
public:
	InstanceNode(void){node_type = INSTANCE_NODE; id = "empty"; type= "empty";};
	char*	id;
	char*	type;
	int		node_type;
	
	OwnSlotsNode own_slots;
	void addOwnSlots(OwnSlotsNode* own_slots);
	virtual void printNode();
	virtual int qualifyNode();
};

void InstanceNode::addOwnSlots(OwnSlotsNode* sn)
{
	own_slots.addOwnSlot(&sn->own_slot);
}



void InstanceNode::printNode()
{
	cout<<"InstanceNode: Node node_type is "<<node_type<<"\n";
	cout<<"InstanceNode: Node type is "<<type<<"\n";
	cout<<"InstanceNode: Node id is "<<id<<"\n";
	own_slots.printNode();
}
int InstanceNode::qualifyNode()
{return node_type;}

class InstancesNode: public Node
{
public:
	InstancesNode(void){node_type = INSTANCES_NODE;};
	int		node_type;
	NodeGroup instances;

	void addInstance(InstanceNode* in);
	virtual void printNode();
	virtual int qualifyNode();
};

void InstancesNode::addInstance(InstanceNode* in)
{
	instances.Sons.push_back(new InstanceNode);
	InstanceNode* ilp = (InstanceNode*) instances.Sons.front();

	ilp->addOwnSlots(&in->own_slots);

}

void InstancesNode::printNode()
{
	cout<<"InstancesNode: Node type is "<<node_type<<"\n";
	list<Node*>::iterator it;
	for(it = instances.Sons.begin(); it!=instances.Sons.end(); ++it)
		(*it)->printNode();
}

int InstancesNode::qualifyNode()
{return node_type;}
