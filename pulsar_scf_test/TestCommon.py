import pulsar as psr

def make_wf():
    MyU=psr.AtomSetUniverse()
    angstrom_to_bohr = 1 / 0.52917721092
    Oy=-0.07579*angstrom_to_bohr
    Hx=0.86681*angstrom_to_bohr
    Hy=0.60144*angstrom_to_bohr
    cg = psr.ShellType.CartesianGaussian
    O1s = psr.BasisShellInfo(cg,0,3,1,
                         [130.709320000, 23.808861000, 6.443608300],
                         [0.15432897, 0.53532814, 0.44463454])
    O2s = psr.BasisShellInfo(cg,0,3,1,
                         [5.033151300, 1.169596100, 0.380389000],
                         [-0.09996723, 0.39951283, 0.70011547])
    O2p = psr.BasisShellInfo(cg,1,3,1,
                         [5.033151300, 1.169596100, 0.380389000],
                         [0.15591627, 0.60768372, 0.39195739])
    H1s = psr.BasisShellInfo(cg,0,3,1,
                         [3.425250910, 0.623913730, 0.168855400],
                         [0.15432897, 0.53532814, 0.44463454])
    O=psr.create_atom([0.00000,Oy,0.00000],8)
    BI=psr.BasisInfo()
    BI.shells=[O1s,O2s,O2p]
    O.basis_sets={"PRIMARY":BI,"FITTING":BI}
    H1=psr.create_atom([Hx,Hy,0.00000],1)
    BI=psr.BasisInfo()
    BI.shells=[H1s]
    H1.basis_sets={"PRIMARY":BI,"FITTING":BI}
    H2=psr.create_atom([-Hx,Hy,0.00000],1)
    BI=psr.BasisInfo()
    BI.shells=[H1s]
    H2.basis_sets={"PRIMARY":BI,"FITTING":BI}
    MyU.insert(O)
    MyU.insert(H1)
    MyU.insert(H2)
    wf=psr.Wavefunction()
    wf.system=psr.System(MyU,True)
    return wf

#MyU=psr.AtomSetUniverse()
#wf=psr.Wavefunction()
#tempsys=psr.make_system("""
#   O 0.0 -0.07579 0.0
#   H 0.86681 0.60144 0.0
#   H -0.86681 0.60144 0.0
# """)
#
#wf.system=psr.apply_single_basis_set("PRIMARY","sto-3g",tempsys)
#return wf
