import ROOT
import copy
from array import array

class Node:
    def __init__(self, name='', formula='', leaves=[]):
        self.name    = name
        self.formula = formula
        self.leaves = leaves

    def create(self):
        formula = copy.deepcopy(self.formula)
        ## format formula
        for i,leaf in enumerate(self.leaves):
            tag = '{'+leaf.name+'}'
            if not tag in formula:
                raise StandardError(leaf.name+' is declared as leaf of node '+self.name+' but it is not part of the given formula.')
            formula = formula.replace(tag,'x[{}]'.format(i))
        ## create TFormula
        tformula = ROOT.TFormula(self.name, formula)
        wtformula = ROOT.WrapperTFormula(tformula, self.name)
        return wtformula

    def replace(self, other):
        self.name    = other.name
        self.formula = other.formula
        self.leaves  = other.leaves

    def find(self, name):
        if self.name==name:
            return self
        else:
            for leaf in self.leaves:
                found = leaf.find(name)
                if found: return found
        return None

    def __str__(self):
        return self.name

    def __eq__(self,other):
        return  self.name==other.name


class Leaf:
    def __init__(self, name='', file='', object='', vars=[]):
        self.name = name
        self.file = file
        self.object = object
        self.vars = vars

    def create(self):
        file = ROOT.TFile.Open(self.file)
        if not file:
            raise StandardError('Cannot open file '+self.file+' for leaf '+self.name)
        object = file.Get(self.object)
        if not object:
            raise StandardError('Cannot load object '+self.object+' for leaf '+self.name)
        ## create wrapper according to object type
        wobject = None
        if isinstance(object, ROOT.TGraph):
            wobject = ROOT.WrapperTGraph(object, self.name)
        elif isinstance(object, ROOT.TH1D):
            wobject = ROOT.WrapperTH1D(object, self.name)
        elif isinstance(object, ROOT.TH2F):
            wobject = ROOT.WrapperTH2F(object, self.name)
        elif isinstance(object, ROOT.TH2D):
            wobject = ROOT.WrapperTH2D(object, self.name)
        elif isinstance(object, ROOT.TH3D):
            wobject = ROOT.WrapperTH3D(object, self.name)
        else:
            raise StandardError('Undefined wrapper for object of class '+str(object.__class__))
        file.Close()
        return wobject

    def replace(self, other):
        self.name = other.name
        self.file = other.file
        self.object = other.object
        self.vars = other.vars

    def find(self, name):
        if self.name==name: 
            return self
        else: return None

    def __str__(self):
        return self.name

    def __eq__(self,other):
        return  self.name==other.name

def flatten(node):
    nodes = []
    if isinstance(node, Node):
        for leaf in node.leaves:
            nodes.extend(flatten(leaf))
        nodes.append(node)
    elif isinstance(node, Leaf):
        nodes.append(node)
    return nodes

def find_node(node, name):
    nodes = flatten(node)
    for n in nodes:
        if n.name==name:
            return n
    return None

def print_tree(node, depth=0):
    string = ''
    for d in xrange(depth-1): string += '  '
    if depth>0: string += '| '
    print string
    string = '  '*(depth-1)+'\_' if depth>0 else ''
    string += str(node)
    print string
    if isinstance(node, Node):
        for leaf in node.leaves:
            print_tree(leaf, depth+1)


class FakeFactor:
    def __init__(self, vars):
        self.fakefactor = ROOT.FakeFactor()
        for var in vars:
            self.fakefactor.addInput(var)
        self.vars = vars
        self.wrapperlist = []
        self.nodelist = []
        self.sys_nodelist = {'':[]}

    def addNode(self, node, sys=''):
        wrapper = node.create()
        leaves = [self.sys_nodelist[sys].index(l) for l in node.leaves] if isinstance(node, Node) else []
        vars = [self.vars.index(v) for v in node.vars] if isinstance(node, Leaf) else []
        self.fakefactor.addNode(wrapper,
                                len(leaves), array('L', leaves if len(leaves)>0 else [0]),
                                len(vars), array('L', vars if len(vars)>0 else [0]),
                                sys
                               )
        if node not in self.nodelist:
            self.wrapperlist.append(wrapper)
            self.nodelist.append(node)
        if node not in self.sys_nodelist[sys]:
            self.sys_nodelist[sys].append(node)

    def systematics(self):
        return self.fakefactor.systematics()

    def value(self, inputs, sys=''):
        return self.fakefactor.value(len(inputs), array('d', inputs), sys)



def fill(fakefactor, node, sys=''):
    if sys not in fakefactor.sys_nodelist: 
        fakefactor.sys_nodelist[sys] = []
        if sys!='': fakefactor.fakefactor.registerSystematic(sys)
    if isinstance(node, Node):
        for leaf in node.leaves:
            fill(fakefactor, leaf, sys)
        fakefactor.addNode(node, sys)
    elif isinstance(node, Leaf):
        fakefactor.addNode(node, sys)
    else:
        raise StandardError('Incompatible fake factor node/leaf type '+str(node.__class__))


def replace(node, name, newnode):
    if isinstance(node, Node):
        # if this is the node to be replaced
        if node.name==name:
            node.replace(newnode)
            return
        # else replace name in formula and apply replace on leaves
        node.formula = node.formula.replace('{'+name+'}', '{'+newnode.name+'}')
        for i,leaf in enumerate(node.leaves):
            if leaf.name==name:
                node.leaves[i] = newnode
            else:
                replace(leaf, name, newnode)
    elif isinstance(node, Leaf):
        if node.name==name:
            node.replace(newnode)
            return

def replace_nodes(root, replacements):
    rootcopy = copy.deepcopy(root)
    for name, newnode in replacements.items():
        replace(rootcopy, name, newnode)
    return rootcopy

        
