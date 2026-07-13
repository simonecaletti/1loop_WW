Unprotect[TypeOf];
Remove[TypeOf];
AppendTo[$Path, "/afs/cern.ch/user/x/xiaol/research/mmapkg/external/LiteRed/Setup"];
AppendTo[$Path, "/afs/cern.ch/user/x/xiaol/research/mmapkg/external/finiteflow-mathtools-master/packages"];
AppendTo[$Path, "/afs/cern.ch/user/x/xiaol/research/mmapkg/external/finiteflow"];
AppendTo[$Path, "/Users/andre/Library/Wolfram/Applications/AMFlow/ibp_interface/FiniteFlow+LiteRed/"];
AppendTo[$LibraryPath, "/afs/cern.ch/user/x/xiaol/research/local/ffmath"];
If[$FrontEnd===Null,$InputFileName,NotebookFileName[]]//DirectoryName//SetDirectory;
If[FileExistsQ["trivial"], Put[{{}, {}}, "results/table"]; Abort[]];
If[!DirectoryQ["ibps"], Abort[]];
<<LiteIBP.m;
<<sup.m;
<<topology.wl;


all = LIBPGetAllInts[FileNames["ibps/ints_tt_*.mx"]];
target = Get["target"];
If[target === 0,
target = Select[all, LIBPIntS[#]<=3 && ((LIBPIntR[#]-LIBPIntT[#])<=0)&],
target = Select[LIBPEliminateZeroSectors[target], #=!=0&]];
If[target === {}, Print["target is empty."]; Put[{{},{}}, "results/table"]; Abort[]];


exints = Import["ibps/exints_def.m"];
all = Join[all, First/@exints];
eqsfiles = FileNames["ibps/sids_*.json"];
vars = {eps};
LIBPWriteSystemJSON[eqsfiles, all, target, vars, "FileName"->"system.json"];


FFNewGraph[graph,in,vars];
FFAlgJSONSparseSolver[graph,ibps,{in},"system.json"];
FFSolverOnlyHomogeneous[graph,ibps];
FFGraphOutput[graph,ibps];
ibplearn = FFSparseSolverLearn[graph,all];
{nonmis, mis} = {"DepVars", "IndepVars"}/.ibplearn;
Print["number of master integrals" -> Length[mis]];
mastersexpr = mis/.exints;


If[mis === {}, Put[{{},{}}, "results/table"]; Abort[]];


FFSparseSolverMarkAndSweepEqs[graph,ibps];
FFSparseSolverDeleteUnneededEqs[graph,ibps];
Print["number of independent eqs."->FFSolverNIndepEqs[graph,ibps]];
sol = FFReconstructFunction[graph,vars,"MaxPrimes"->200,"MaxDegree"->1000];
sol = Partition[sol, Length@mis]//Factor;
sol = Join[sol, IdentityMatrix[Length@mis]];
rule = Thread[Join[nonmis, mastersexpr]->sol];
Put[{mastersexpr, rule}, "results/table"];


FFDeleteGraph[graph];
Quit[];