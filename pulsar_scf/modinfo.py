from pulsar import OptionType

#Stuff common to many modules
c_mod,base,modpath="c_module","MatrixFactory","pulsar_scf.so"
version="0.1a"
no_options={}
no_ref,Ryan=[""],["Ryan Richard"]
mods=["TElectronic","NuclearElectronic","HCore","Overlap","JK",
      "DFJK","G","F","SchwarzScreen","Orthogonalizer","Metric",
      "GAJK","GAT","GAV","GAS","GAX","GAH","GAG","GAF"
      ]
force_cache=(OptionType.Bool,False,False,None,
    "Should this module fail if the cache is not used?")
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
    "METRIC_INTS_KEY":(OptionType.String,None,True,None,"The key for the metric integrals"),
    "FORCE_CACHE":force_cache
}
minfo["TElectronic"]["description"]="Builds the kinetic energy in the AO basis"
minfo["TElectronic"]["options"]={
   "T_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
   "FORCE_CACHE":force_cache
   }
minfo["GAT"]["description"]="Builds the kinetic energy in the AO basis"
minfo["GAT"]["options"]={
   "T_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
   "FORCE_CACHE":force_cache
    }
minfo["NuclearElectronic"]["description"]="Builds the nucleus-electron energy in the AO basis"
minfo["NuclearElectronic"]["options"]={
      "V_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
      "FORCE_CACHE":force_cache
   }
minfo["GAV"]["description"]="Builds the nucleus-electron energy in the AO basis"
minfo["GAV"]["options"]={
      "V_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
      "FORCE_CACHE":force_cache
 }
minfo["HCore"]["description"]="Builds the core Hamiltonian in the AO basis"
minfo["HCore"]["options"]={
      "H_KEYS":(OptionType.ListString,["PSR_T","PSR_V"],False,None,"The keys for the core Hamiltonian"),
      "FORCE_CACHE":force_cache
   }
minfo["GAH"]["description"]="Builds the core Hamiltonian in the AO basis"
minfo["GAH"]["options"]={
      "H_KEYS":(OptionType.ListString,["PSR_GAT","PSR_GAV"],False,None,"The keys for the core Hamiltonian"),
      "FORCE_CACHE":force_cache
    }
minfo["Overlap"]["description"]="Builds the overlap matrix in the AO basis"
minfo["Overlap"]["options"]={
      "S_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
      "FORCE_CACHE":force_cache
   }
minfo["GAS"]["description"]="Builds the overlap matrix in the AO basis"
minfo["GAS"]["options"]={
         "S_INTS_KEY":(OptionType.String,None,True,None,"The key for the integrals module"),
         "FORCE_CACHE":force_cache
      }
minfo["JK"]["description"]="Builds the J and K matrices"
minfo["JK"]["options"]={
      "ERI_KEY":(OptionType.String,None,True,None,"The key for the ERI"),
      "FORCE_CACHE":force_cache
   }
minfo["DFJK"]["description"]="Builds the J and K matrices using density fitting"
minfo["DFJK"]["options"]={
    "FITTING_BASIS_KEY":(OptionType.String,"FITTING",False,None,
        "The key for the fitting basis set"),
    "FITTING_COEF_KEY":(OptionType.String,"PSR_DFCoef",False,None,
        "The key for module that builds the density fitting coeficients"),
    "FORCE_CACHE":force_cache
}
minfo["GAJK"]["description"]="Builds the J and K matrices in a direct manner using Global Arrays"
minfo["GAJK"]["options"]={
    "ERI_KEY":(OptionType.String,None,True,None,"The key for the ERI"),
    "FORCE_CACHE":force_cache
}
minfo["G"]["description"]="Builds G"
minfo["G"]["options"]={
         "JK_KEY":(OptionType.String,"PSR_JK",False,None,"The key for the JK builder"),
         "FORCE_CACHE":force_cache
      }
minfo["GAG"]["description"]="Builds G"
minfo["GAG"]["options"]={
          "JK_KEY":(OptionType.String,"PSR_GA_JK",False,None,"The key for the JK builder"),
          "FORCE_CACHE":force_cache
     }
minfo["F"]["description"]="Builds the Fock Matrix"
minfo["F"]["options"]={
               "G_KEY":(OptionType.String,"PSR_G",False,None,"The key for the G builder"),
               "H_KEY":(OptionType.String,"PSR_H",False,None,"The key for the H builder"),
               "FORCE_CACHE":force_cache
            }
minfo["GAF"]["description"]="Builds the Fock Matrix"
minfo["GAF"]["options"]={
                "G_KEY":(OptionType.String,"PSR_GAG",False,None,"The key for the G builder"),
                "H_KEY":(OptionType.String,"PSR_GAH",False,None,"The key for the H builder"),
                "FORCE_CACHE":force_cache
            }
minfo["SchwarzScreen"]["description"]="Builds diagonal elements of ERI"
minfo["SchwarzScreen"]["options"]={
    "ERI_INTS_KEY":(OptionType.String,"PSR_ERI",False,None,"The key for ERI builder."),
    "FORCE_CACHE":force_cache
}
minfo["Orthogonalizer"]["description"]="Builds a transformation that orthogonalizes the AOs"
minfo["Orthogonalizer"]["options"]={
    "S_KEY":(OptionType.String,"PSR_S",False,None,"The key for building the overlap matrix"),
    "FORCE_CACHE":force_cache
}
minfo["GAX"]["description"]="Builds a transformation that orthogonalizes the AOs"
minfo["GAX"]["options"]={
    "S_KEY":(OptionType.String,"PSR_GAS",False,None,"The key for building the overlap matrix"),
    "FORCE_CACHE":force_cache
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
    "FORCE_CACHE":force_cache
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
    "X_KEY":(OptionType.String,"PSR_X",False,None,"The key for the X builder"),
    "FORCE_CACHE":force_cache
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
    "FORCE_CACHE":force_cache
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
        "The key for three-center, two-electron integral tensor builder"),
    "FORCE_CACHE":force_cache
  }
}
