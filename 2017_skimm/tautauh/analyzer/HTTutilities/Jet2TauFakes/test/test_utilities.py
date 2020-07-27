from HTTutilities.Jet2TauFakes.Utilities import Leaf, Node, fill, FakeFactor,replace_nodes 
import ROOT
from array import array


print '########################################'
print 'Test 1'
ff = FakeFactor(vars=['tau_pt', 'tau_decay', 'mt'])

leaf1 = Leaf(
    name='ff_W',
    file='/afs/cern.ch/work/j/jsauvan/public/HTauTau/FakeFactors/FakeFactors_Data_HighMT_2D.root',
    object='FakeFactors_Data_HighMT_2D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode',
    vars=['tau_pt','tau_decay']
)
leaf2 = Leaf(
    name='corr_mt',
    file='/afs/cern.ch/work/j/jsauvan/public/HTauTau/FakeFactors/mtCorrections.root',
    object='mt_correction',
    vars=['mt']
)
leaf2_sys = Leaf(
    name='corr_mt_statup',
    file='/afs/cern.ch/work/j/jsauvan/public/HTauTau/FakeFactors/mtCorrections.root',
    object='mt_correction_statup',
    vars=['mt']
)
node = Node(
    name='ff_W_corr',
    formula='{corr_mt}*{ff_W}',
    leaves=[leaf1,leaf2]
)
node_sys = Node(
    name='ff_W_corr_statup',
    formula='({corr_mt} + {corr_mt_statup} - {corr_mt})*{ff_W}',
    leaves=[leaf1,leaf2,leaf2_sys]
)

fill(ff, node)
nodesys = replace_nodes(node, {'corr_mt':leaf2_sys})
nodesys2 = replace_nodes(node, {'ff_W_corr':node_sys})
fill(ff, nodesys, 'Sys_CorrMT_StatUp')
fill(ff, nodesys2, 'Sys_CorrMT_StatUp2')
ff_sys = FakeFactor(vars=['tau_pt', 'tau_decay', 'mt'])
ff_sys2 = FakeFactor(vars=['tau_pt', 'tau_decay', 'mt'])
fill(ff_sys, nodesys)
fill(ff_sys2, nodesys2)

print ff.value([30,1,10])
print ''
print ff_sys.value([30,1,10])
print ff_sys2.value([30,1,10])
print ''
print ff.value([30,1,10], 'Sys_CorrMT_StatUp')
print ff.value([30,1,10], 'Sys_CorrMT_StatUp2')


print '########################################'
print 'Test 2'
ff2 = FakeFactor(vars=['tau_pt', 'tau_decay', 'mt'])
root = Node(
    name='ff_W_corr',
    formula='{corr_mt}*{ff_W}',
    leaves=[
        Leaf(
            name='ff_W',
            file='/afs/cern.ch/work/j/jsauvan/public/HTauTau/FakeFactors/FakeFactors_Data_HighMT_2D.root',
            object='FakeFactors_Data_HighMT_2D_Iso_Medium_InvertIso_Medium_tau_pt_vs_decayMode',
            vars=['tau_pt','tau_decay']
        ),
        Leaf(
            name='corr_mt',
            file='/afs/cern.ch/work/j/jsauvan/public/HTauTau/FakeFactors/mtCorrections.root',
            object='mt_correction',
            vars=['mt']
        )
    ]
)
fill(ff2, root)
print ff2.value([30,1,10])
