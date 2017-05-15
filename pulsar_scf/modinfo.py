from pulsar import OptionType

#Stuff common to many modules
c_mod,base,modpath="c_module","MatrixFactory","pulsar_scf.so"
version="0.1a"
no_options={}
no_ref,Ryan=[""],["Ryan Richard"]
mods=["TElectronic","NuclearElectronic","HCore","Overlap","JK",
      "DFJK","G","F","SchwarzScreen","Orthogonalizer","Metric"]
minfo={}
for mod in mods:
    minfo[mod]={
            "type":c_mod,
            "base":base,
            "modpath":modpath,
            "version":version,
            "authors":Ryan,
            "refs":no_ref,
            "options":no_options}
minfo["Metric"]["description"]="Builds the density fitting metric matrix"
minfo["Metric"]["options"]={
    "METRIC_INTS_KEY":(OptionType.String,None,True,None,"The key for the metric integrals")
}
minfo["TElectronic"]["description"]="Builds the kinetic energy in the AO basis"
minfo["TElectronic"]["options"]={
   "T_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
   }
minfo["NuclearElectronic"]["description"]="Builds the nucleus-electron energy in the AO basis"
minfo["NuclearElectronic"]["options"]={
      "V_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
   }
minfo["HCore"]["description"]="Builds the core Hamiltonian in the AO basis"
minfo["HCore"]["options"]={
      "H_KEYS":(OptionType.ListString,["PSR_T","PSR_V"],False,None,"The keys for the core Hamiltonian"),
   }
minfo["Overlap"]["description"]="Builds the overlap matrix in the AO basis"
minfo["Overlap"]["options"]={
      "S_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
   }
minfo["JK"]["description"]="Builds the J and K matrices"
minfo["JK"]["options"]={
      "ERI_KEY":(OptionType.String,None,True,None,"The key for the ERI"),
   }
minfo["DFJK"]["description"]="Builds the J and K matrices using density fitting"
minfo["DFJK"]["options"]={
    "FITTING_BASIS_KEY":(OptionType.String,"FITTING",False,None,
        "The key for the fitting basis set"),
    "FITTING_COEF_KEY":(OptionType.String,"PSR_DFCoef",False,None,
        "The key for module that builds the density fitting coeficients")
}
minfo["G"]["description"]="Builds G"
minfo["G"]["options"]={
         "JK_KEY":(OptionType.String,"PSR_JK",False,None,"The key for the JK builder"),
      }
minfo["F"]["description"]="Builds the Fock Matrix"
minfo["F"]["options"]={
               "G_KEY":(OptionType.String,"PSR_G",False,None,"The key for the G builder"),
               "H_KEY":(OptionType.String,"PSR_H",False,None,"The key for the H builder")
            }
minfo["SchwarzScreen"]["description"]="Builds diagonal elements of ERI"
minfo["SchwarzScreen"]["options"]={
    "ERI_INTS_KEY":(OptionType.String,"PSR_ERI",False,None,"The key for ERI builder.")
}
minfo["Orthogonalizer"]["description"]="Builds a transformation that orthogonalizes the AOs"
minfo["Orthogonalizer"]["options"]={
    "S_KEY":(OptionType.String,"PSR_S",False,None,"The key for building the overlap matrix")
}
minfo["CoreDensity"]={
  "type":c_mod,
  "base":"EnergyMethod",
  "modpath":modpath,
  "version":version,
  "authors":Ryan,
  "refs":no_ref,
  "description":"Returns the energy and wavefunction for a set of non-interacting electrons",
  "options":{
    "H_KEY":(OptionType.String,"PSR_H",False,None,"The key for the H builder"),
    "S_KEY":(OptionType.String,"PSR_S",False,None,"The key for the S builder"),
  }
}


minfo["SCF"]={
  "type":c_mod,
  "base":"EnergyMethod",
  "modpath":modpath,
  "version":version,
  "authors":Ryan,
  "refs":no_ref,
  "description":"Runs an SCF computation",
  "options":{
    "H_KEY":(OptionType.String,"PSR_H",False,None,"The key for the H builder"),
    "F_KEY":(OptionType.String,"PSR_F",False,None,"The key for the F builder"),
    "S_KEY":(OptionType.String,"PSR_S",False,None,"The key for the S builder"),
    "X_KEY":(OptionType.String,"PSR_X",False,None,"The key for the X builder")
  }
}


minfo["DFInts"]={
  "type":c_mod,
  "base":"Rank3Builder",
  "modpath":modpath,
  "version":version,
  "authors":Ryan,
  "refs":no_ref,
  "description":"Builds the 3 center, 2 e integrals",
  "options":{
    "DF_INTS_KEY":(OptionType.String,None,True,None,"The key for the DF integals"),
  }
}
minfo["DFCoef"]={
  "type":c_mod,
  "base":"Rank3Builder",
  "modpath":modpath,
  "version":version,
  "authors":Ryan,
  "refs":no_ref,
  "description":"Builds the 3 center, 2 e integrals",
  "options":{
    "METRIC_KEY":(OptionType.String,"PSR_Metric",False,None,
        "The key for metric matrix builder"),
    "DF_INTS_KEY":(OptionType.String,"PSR_3C2E",False,None,
        "The key for three-center, two-electron integral tensor builder")
  }
}
