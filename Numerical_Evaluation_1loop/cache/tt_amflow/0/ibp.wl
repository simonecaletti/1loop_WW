Unprotect[TypeOf];
Remove[TypeOf];
AppendTo[$Path, "/afs/cern.ch/user/x/xiaol/research/mmapkg/external/LiteRed/Setup"];
AppendTo[$Path, "/afs/cern.ch/user/x/xiaol/research/mmapkg/external/finiteflow-mathtools-master/packages"];
AppendTo[$Path, "/afs/cern.ch/user/x/xiaol/research/mmapkg/external/finiteflow"];
AppendTo[$Path, "/Users/andre/Library/Wolfram/Applications/AMFlow/ibp_interface/FiniteFlow+LiteRed/"];
AppendTo[$LibraryPath, "/afs/cern.ch/user/x/xiaol/research/local/ffmath"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
<<LiteIBP.m;
<<sup.m;
<<topology.wl;


If[NonZeroSectors[tt]==={}, Print["this family is trivial."]; Put[Null, "trivial"]; Quit[]];


rank=3;
maxdots=0;
GetSeeds["IBP"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots},0,rank];
GetSeeds["LI"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots},0,rank];
GetSeeds["SR"][sector_]:=LIBPGenerateSeeds[sector,{0,0},0,rank];
GetSeeds["Map"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots+1},0,rank];
GetSeeds["ExtMap"][sector_]:=LIBPGenerateSeeds[sector,{0,maxdots+1},0,rank];
LIBPFastGenIds[tt,GetSeeds,"Directory"->"ibps","LaunchKernels"->4];


exints = {};
Export["ibps/exints_def.m",exints];
Export["ibps/ids_tt_exints.mx",(-#[[1]]+#[[2]])&/@exints,"MX"];
Export["ibps/ints_tt_exints.mx",LIBPIntegralsIn[exints],"MX"];


LIBPSerializeFastIds[FileNames["ibps/ids_tt_*.mx"], {s -> 30, t -> -10/3, msq -> 1}, {eps}, "ExtraInts"->(First/@exints)];


Quit[];