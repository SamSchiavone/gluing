/*
  AttachSpec("~/github/quartic_reconstruction/magma/spec");
  AttachSpec("~/github/curve_reconstruction/magma/spec");
  AttachSpec("~/github/endomorphisms/endomorphisms/magma/spec");
  AttachSpec("~/github/quartic_isomorphisms/magma/spec");
*/

AttachSpec("~/github/CHIMP/quartic_reconstruction/magma/spec");
AttachSpec("~/github/CHIMP/curve_reconstruction/magma/spec");
AttachSpec("~/github/CHIMP/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/github/CHIMP/quartic_isomorphisms/magma/spec");
AttachSpec("spec");

/* Arithmetic reconstruction */

SetVerbose("QuarticIso", 1);
SetVerbose("QuarticRec", 1);
SetVerbose("Gluing", 1);
//SetVerbose("CurveRec", 1);
SetVerbose("EndoFind", 1);

//prec := 500;
prec := 20;
F := RationalsExtra(prec);
R<x> := PolynomialRing(F);

// TODO: breaks for odd degree models
/*
f1 := x^6 - x;
f2 := x^6 + 20*x^3 + 36*x;
*/

f1 := &*[x-el : el in [0..5]];
f2 := &*[x-el : el in [10..15]];

print "";
print "Can we glue over QQ?";
print ExistsGaloisStableSubgroupFor22(f1, f2);

X1 := HyperellipticCurve(f1);
X2 := HyperellipticCurve(f2);

print "";
print "All arithmetic 2-gluings:";
Ys := AllArithmetic2GluingsCCFor22(X1, X2, F);

Pi := Ys[1];
tau := SmallPeriodMatrix(Pi);
// now need to reduce...see reduce.m line 258

// TODO: adapt for genus 4
/* Reduce small period matrix */
//taunew, gamma := ReduceSmallPeriodMatrix(tau);
//Imtaunew := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(taunew) ]);
//
//vprint CurveRec, 2 : "";
//vprint CurveRec, 2 : "Eigenvalues of imaginary part of reduced tau:";
//vprint CurveRec, 2 : [ ComplexField(5) ! tup[1] : tup in Eigenvalues(Imtaunew) ];
//
///* Calculate corresponding big period matrix */
//A := Submatrix(gamma, 1,1, 2,2);
//B := Submatrix(gamma, 1,3, 2,2);
//C := Submatrix(gamma, 3,1, 2,2);
//D := Submatrix(gamma, 3,3, 2,2);
//Pnew := P * Transpose(BlockMatrix([[A, B], [C, D]]));
//P1new := Submatrix(Pnew, 1,1, 2,2);
//P2new := Submatrix(Pnew, 1,3, 2,2);
//P2inew := P2new^(-1);

/* Calculation of theta derivatives at odd two-torsion points */
char_pairs := [];
for epsilon0, eta0 in CartesianPower([0,1],4) do
  epsilon := Matrix(GF(2),1,4,[el : el in epsilon0]);
  eta := Matrix(GF(2),1,4,[el : el in eta0]);
  if (epsilon*Transpose(eta))[1,1] eq 1 then // dot product = 1
    Append(~char_pairs, [epsilon, eta]);
  end if;
end for;
assert #char_pairs eq 120;

taunew := tau;
CC := BaseRing(Parent(taunew));
ws := [[* pair, (1/2)*taunew*Transpose(Matrix(CC,1,4, ChangeUniverse(Eltseq(pair[1]),ZZ))) + (1/2)*Transpose(Matrix(CC,1,4, ChangeUniverse(Eltseq(pair[2]),ZZ))) *] : pair in char_pairs];

vprint CurveRec : "";
vprint CurveRec : "Calculating theta derivatives...";
//theta_derss := [ ThetaDerivatives(taunew, w[2]) : w in ws ];
theta_derss := [ ThetaDerivatives(taunew, w[2] : B := 5) : w in ws ];
vprint CurveRec : "done calculating theta derivatives.";

function ArfInvariant(pair)
  return (pair[1]*Transpose(pair[2]))[1,1];
end function;

function IsSyzygetic(pair1,pair2)
  inv := GF(2)!0;
  for el in pair1 do
    inv +:= ArfInvariant(el);
  end for;
  for el in pair1 do
    inv +:= ArfInvariant(el);
  end for;
  return (inv eq 0);
end function;

function IsSteinerCompatible(pair1,pair2)
  if #(pair1 meet pair2) ne 0 then
    return false;
  end if;
  return IsSyzygetic(pair1,pair2);
end function;

print "Computing Steiner systems";
theta_set := SequenceToSet(char_pairs);
pairs := Subsets(theta_set,2);
blocks := [[] : i in [1..255]];
for pair in pairs do
  for block in blocks do
    found_bool := true;
    for el in block do
      if not IsSteinerCompatible(pair,el) then
        found_bool := false;
        break;
      end if;
      if found_bool then
        Append(~block, pair);
      else
        Append(~blocks, [pair]);
      end if;
    end for;
  end for;
end for;



/*
  blocks := [];
  for pair in pairs do
    for block in blocks do
      found_bool := true;
      for el in block do
        if not IsSteinerCompatible(pair,el) then
          found_bool := false;
          break;
        end if;
        if found_bool then
          Append(~block, pair);
        else
          Append(~blocks, [pair]);
        end if;
      end for;
    end for;
  end for;
*/


tetrads := [];
for pair1, pair2, pair3 in char_pairs do
  s := [pair1[i] + pair2[i] + pair3[i] : i in [1..#pair1]];
  if ArfInvariant(pair1) + ArfInvariant(pair2) + ArfInvariant(pair3) + ArfInvariant(s) eq 0 then
    Append(~tetrads, [pair1, pair2, pair3, s]);
  end if;
end for;

/* Determination of ratios = roots */
Hs := [ Matrix(CC, [ theta_ders ]) * P2inew : theta_ders in theta_derss ];

S<x,y,z,w> := PolynomialRing(CC,4);
planes := [];
for el in theta_derss do
  Append(~planes, &+[el[i]*S.i : i in [1..4]]);
end for;

/*
  [0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 + 1.58987247748376043097412350231704567737561104926950398528228906437188683644838494217146732983486999852993278786867807625209465536905933249358136967124145403586217417818823765369154365377133378218914649458628866878982525801881675978599117337095386958553540690401340946073930400760727704201639387254245523040650921946536803675644452232643076626970181596840034314074532990174596902773345360706331610218411555133903915130783936480518158446763766439911977681808124206415225977261713555330426529290614409948347573274E-526*I 0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 + 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I 0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 + 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 - 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I 0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 + 1.58987247748376043097412350231704567737561104926950398528228906437188683644838494217146732983486999852993278786867807625209465536905933249358136967124145403586217417818823765369154365377133378218914649458628866878982525801881675978599117337095386958553540690401340946073930400760727704201639387254245523040650921946536803675644452232643076626970181596840034314074532990174596902773345360706331610218411555133903915130783936480518158446763766439911977681808124206415225977261713555330426529290614409948347573274E-526*I 0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 - 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I]
  [0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 1.27934012647429327816867395577872541993952830377772623117819688530861544563372368539856282525096751908426671401312463196667888220198729500902145105553484734456388050702835648209405764931018905900505481341390879272586455608606953553434182115908481763879864386866742681954061468486375470766251372703007864198829316925721977847673433177857807653308085567362103581454840485108956172190103397034200929410367189253970990198236790832649718567164025130167607703449047876251647160774933740600558673199806628213002742626 + 5.49485318457877991359432924774149605119608208964535907118701409819941609493248784770078844107520073453008165611254359023858114739168618837162507094838913430255922399917423651858435992347096891069519003005025893688959882603448641909701089052339687546232505580995204105137644291560044385804533924216328942962670019629832634123479462713573496355756436949435209714082935638614402936438489203184590538274968390991856882996553546990597289030871877363139630057684509585590884603988519158044250857283897435644077064193E-526*I 0.571471522261315178722590973066918653703230213290226828176959768479325552114815546595532239139475367073564804814281912761337967453541270056014731437430893933521795869234004770831726194375887820405106954067428358592307828547448410554938527384997544920088239192360391820780271564544117639586566314942475201271004311447517117843147205754577505524593547881548308801512411639413908464171771396854569191498553813216396666939408425490875321461720536795206778644062512790111747560440280584533810159809549275461829524737 + 0.958886356306729701282398146416919704449977463796235895597114756189288523471117300952413103635019129270326872099839712933645754606252284293310488703013388125645562997286195242898725587671468705502474449845396398726795062475036524599845592550067209197654126562967029857851409333778940141836825309174099374287426735413551859670917399976870845816894722334745760181808600788526340639059934897945849748015935984081819677347538092397631754016593690499490210841951305590477961176923051954964236540712352899259291445916*I 0.571471522261315178722590973066918653703230213290226828176959768479325552114815546595532239139475367073564804814281912761337967453541270056014731437430893933521795869234004770831726194375887820405106954067428358592307828547448410554938527384997544920088239192360391820780271564544117639586566314942475201271004311447517117843147205754577505524593547881548308801512411639413908464171771396854569191498553813216396666939408425490875321461720536795206778644062512790111747560440280584533810159809549275461829524737 + 0.130567271188432618058163959950389209447966651798401543265052224231509770709131324988249391487577695775438057582290468015348975440397891168734290300608281301423756822543799057064147605560605480206960200759483916019903744423363642216028965881271331045508409624658528089123919985315912014175015034359681811578304709802693607294569390573407111075408992572079024787892162605715226981670349151068764438603862807892688120283052737576762808378562641842723897837264415253085100335532190579640004489449754466659614934633*I 0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 0.571471522261315178722590973066918653703230213290226828176959768479325552114815546595532239139475367073564804814281912761337967453541270056014731437430893933521795869234004770831726194375887820405106954067428358592307828547448410554938527384997544920088239192360391820780271564544117639586566314942475201271004311447517117843147205754577505524593547881548308801512411639413908464171771396854569191498553813216396666939408425490875321461720536795206778644062512790111747560440280584533810159809549275461829524737 - 0.958886356306729701282398146416919704449977463796235895597114756189288523471117300952413103635019129270326872099839712933645754606252284293310488703013388125645562997286195242898725587671468705502474449845396398726795062475036524599845592550067209197654126562967029857851409333778940141836825309174099374287426735413551859670917399976870845816894722334745760181808600788526340639059934897945849748015935984081819677347538092397631754016593690499490210841951305590477961176923051954964236540712352899259291445916*I 1.27934012647429327816867395577872541993952830377772623117819688530861544563372368539856282525096751908426671401312463196667888220198729500902145105553484734456388050702835648209405764931018905900505481341390879272586455608606953553434182115908481763879864386866742681954061468486375470766251372703007864198829316925721977847673433177857807653308085567362103581454840485108956172190103397034200929410367189253970990198236790832649718567164025130167607703449047876251647160774933740600558673199806628213002742626 + 5.49485318457877991359432924774149605119608208964535907118701409819941609493248784770078844107520073453008165611254359023858114739168618837162507094838913430255922399917423651858435992347096891069519003005025893688959882603448641909701089052339687546232505580995204105137644291560044385804533924216328942962670019629832634123479462713573496355756436949435209714082935638614402936438489203184590538274968390991856882996553546990597289030871877363139630057684509585590884603988519158044250857283897435644077064193E-526*I 0.571471522261315178722590973066918653703230213290226828176959768479325552114815546595532239139475367073564804814281912761337967453541270056014731437430893933521795869234004770831726194375887820405106954067428358592307828547448410554938527384997544920088239192360391820780271564544117639586566314942475201271004311447517117843147205754577505524593547881548308801512411639413908464171771396854569191498553813216396666939408425490875321461720536795206778644062512790111747560440280584533810159809549275461829524737 - 0.130567271188432618058163959950389209447966651798401543265052224231509770709131324988249391487577695775438057582290468015348975440397891168734290300608281301423756822543799057064147605560605480206960200759483916019903744423363642216028965881271331045508409624658528089123919985315912014175015034359681811578304709802693607294569390573407111075408992572079024787892162605715226981670349151068764438603862807892688120283052737576762808378562641842723897837264415253085100335532190579640004489449754466659614934633*I]
  [0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 + 1.58987247748376043097412350231704567737561104926950398528228906437188683644838494217146732983486999852993278786867807625209465536905933249358136967124145403586217417818823765369154365377133378218914649458628866878982525801881675978599117337095386958553540690401340946073930400760727704201639387254245523040650921946536803675644452232643076626970181596840034314074532990174596902773345360706331610218411555133903915130783936480518158446763766439911977681808124206415225977261713555330426529290614409948347573274E-526*I 3.74342216623218556215132021904732687935330866497899777650060486150503591613325117810646710280597943962261826197083935160892491770783004146193284463151219686876197585284729449422107721080011806576705562527396359010900602621944080953870373079788630301233229570162376363602187290234522519859974022951496383612551645643177415260413765843633853045125111822372839676019677597572909511362400723151642428999133687035118143968375843598598216471204477561220761613194022213120632726133536972322979972911402291535988330274E-526 - 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I -0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 - 1.58987247748376043097412350231704567737561104926950398528228906437188683644838494217146732983486999852993278786867807625209465536905933249358136967124145403586217417818823765369154365377133378218914649458628866878982525801881675978599117337095386958553540690401340946073930400760727704201639387254245523040650921946536803675644452232643076626970181596840034314074532990174596902773345360706331610218411555133903915130783936480518158446763766439911977681808124206415225977261713555330426529290614409948347573274E-526*I -3.74342216623218556215132021904732687935330866497899777650060486150503591613325117810646710280597943962261826197083935160892491770783004146193284463151219686876197585284729449422107721080011806576705562527396359010900602621944080953870373079788630301233229570162376363602187290234522519859974022951496383612551645643177415260413765843633853045125111822372839676019677597572909511362400723151642428999133687035118143968375843598598216471204477561220761613194022213120632726133536972322979972911402291535988330274E-526 + 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I 7.48684433246437112430264043809465375870661732995799555300120972301007183226650235621293420561195887924523652394167870321784983541566008292386568926302439373752395170569458898844215442160023613153411125054792718021801205243888161907740746159577260602466459140324752727204374580469045039719948045902992767225103291286354830520827531687267706090250223644745679352039355195145819022724801446303284857998267374070236287936751687197196432942408955122441523226388044426241265452267073944645959945822804583071976660548E-526 - 0.435781450998064927736224842546923565559177646237854975544866792168319317672099450376264998049038730018305971872852072379597892018660070184817911601448667770827727927931997719985149277292829674283773860241952125898679522759360066726349823372535416097265014475050223178790131727637940862404736137413512474346292578086498186786194716220111182756921485962729913987880305357696627048292113619605845674647919516789803119052236331989757824958062532936885643471686288337425224604982097013841696412064842946367562552219*I -0.370162329747121691378252985769128814728551703413590611871031330757588199549707846398819012878088577231566303765481308945603369931105713013007236498593148255617135275252472250585156768737215375882032353496267430263634476926703589217856069708816472511777376612205563728064177249881574469449816008394510768651859496140947379263976307506631116411534880711033868923212163298100694037214561073439315697120445141151221313784355266763474501426672157619376571135710598310525643833637923598107879378361523111518371390199 - 1.58987247748376043097412350231704567737561104926950398528228906437188683644838494217146732983486999852993278786867807625209465536905933249358136967124145403586217417818823765369154365377133378218914649458628866878982525801881675978599117337095386958553540690401340946073930400760727704201639387254245523040650921946536803675644452232643076626970181596840034314074532990174596902773345360706331610218411555133903915130783936480518158446763766439911977681808124206415225977261713555330426529290614409948347573274E-526*I -3.74342216623218556215132021904732687935330866497899777650060486150503591613325117810646710280597943962261826197083935160892491770783004146193284463151219686876197585284729449422107721080011806576705562527396359010900602621944080953870373079788630301233229570162376363602187290234522519859974022951496383612551645643177415260413765843633853045125111822372839676019677597572909511362400723151642428999133687035118143968375843598598216471204477561220761613194022213120632726133536972322979972911402291535988330274E-526 + 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I 3.74342216623218556215132021904732687935330866497899777650060486150503591613325117810646710280597943962261826197083935160892491770783004146193284463151219686876197585284729449422107721080011806576705562527396359010900602621944080953870373079788630301233229570162376363602187290234522519859974022951496383612551645643177415260413765843633853045125111822372839676019677597572909511362400723151642428999133687035118143968375843598598216471204477561220761613194022213120632726133536972322979972911402291535988330274E-526 - 0.217890725499032463868112421273461782779588823118927487772433396084159658836049725188132499024519365009152985936426036189798946009330035092408955800724333885413863963965998859992574638646414837141886930120976062949339761379680033363174911686267708048632507237525111589395065863818970431202368068706756237173146289043249093393097358110055591378460742981364956993940152678848313524146056809802922837323959758394901559526118165994878912479031266468442821735843144168712612302491048506920848206032421473183781276110*I]
  [4.27309481973253209250512083075820680098874724742613294688727307605520754761189401058372236792036113938922784246909500221737166676459840018608709642336237648969314862175872727668329388174804157922543048903010266122865259781448430273349922447316227003786200531441602910142204406335986233408472639888758288778959927285699091048291028082088866963994235499188699803363404462042084883631738213124772616270300522472860980478296109312562033572844211298897249000116849589536818589681951656561260394342478039064554342672 + 1.83532339776776265510175232870141891094380836493972007080473218673788864517932862848712215064478492434189101619176275535428386831442037990320954741202326762353733886736493282867887948830990378221420073887440710949577800442477649776928567100409111683207060477502491008597131171985087117722005688459735390264709380956821942101118732077607265304879628700316530499703505921400455060516169120741038361500120493601670184396520550021717508767058461923554003516943735351063750505493718417191214708207531871358113653491E-525*I 5.39081511188198310179590041653676667757150727264901776126586031049692562073779132715429917919427486614565251640918644362554342731912273511034790337480920545157575959751419407189957485881553540287409919517949198358871572459569827390042106968164942134433471092167805247092866858446950012358560169453254324242951479205172889272087451533604179717310455722567599516000157078206047638391319142244879511851881128877697801484429497004040128636708210057619271237107355019323892075881554403597788751760918275755820239917E-525 - 3.13779361129705433996352235915153753224586569498551077332144871703088511183161455283373809388021277936185673146410007483163521469955263521740004671025672697978420263694618384282447197413561707692134375105515702822019267627183685823159470941274428968397919893821814575180206797196864445386050599624166174601888962584604279360189098107742675960150215214839533012121012757700947588052050299597507812125553356803083527260871975234642087880690635518391842820038274727760408420183353702417271860103656763109710420701*I -4.98096342394551019195120381347001356722504533791363234988851019288449744113080214938675295403185329139992975166793772142271258151304442513909381604146632990073523325955307898794562533668234281782537834837658309536220932535310542771290251824724954275657240999072306410018238718367949940216067381097518632850688813066669357111649740684488924064842966278395972504667003783209650209404664470473516626530812330405192303982592057596124219993836182749544178839159646186777290994412857338708438051561329739731374132824 - 2.13935779594163842233355642709119528249521925823403989240099047419182844594163372694154617394239007198294095347993243527595277010822795133074387676608036746611809657810566130554997964611843067325866549759131456247878514062226540169569226242329355713176791248500861356587694829916732142782092779675878417336917923909520067087992398504000426262545825291783555285482826554036037196417194281024790664045908394233089603430439291179577887349850954176225940687576575164974314437660565471134851615019004153512755279693E-525*I -5.39081511188198310179590041653676667757150727264901776126586031049692562073779132715429917919427486614565251640918644362554342731912273511034790337480920545157575959751419407189957485881553540287409919517949198358871572459569827390042106968164942134433471092167805247092866858446950012358560169453254324242951479205172889272087451533604179717310455722567599516000157078206047638391319142244879511851881128877697801484429497004040128636708210057619271237107355019323892075881554403597788751760918275755820239917E-525 + 3.13779361129705433996352235915153753224586569498551077332144871703088511183161455283373809388021277936185673146410007483163521469955263521740004671025672697978420263694618384282447197413561707692134375105515702822019267627183685823159470941274428968397919893821814575180206797196864445386050599624166174601888962584604279360189098107742675960150215214839533012121012757700947588052050299597507812125553356803083527260871975234642087880690635518391842820038274727760408420183353702417271860103656763109710420701*I 7.93548060739696160716480026216310104162352877959695335997130368653132833919067323622373715564134746582178627703582387079353773390090473708896841640794257344065836006920808432730623633636951952308707973601083398336759868190580749989267651462613267237299205666476271323825202734278712574582749775850973269576855269805541297757893926150960905791004647666728999348098073831452452280029365331268453121291906177420195116873020223984910825082605967690865265591755401026955379478904576054419322361035174906168301171536E-525 - 4.61894905235751451347857634537001407448770976597535284197877237014621271813925715373914876346554269173393583389310165982667687106739648418564769661570324031112479292440757531397978798404950770325165900393848909102660271644032795169555616548789682306366696399981928796614915724701123265239739144285448836661953520047036908245108594334792604972003284477145718945458737878839672444626183449819598562368692078368340743108846879505110386633775061305430423039139171388042244672088535129769697309954793839699485539146*I -4.98096342394551019195120381347001356722504533791363234988851019288449744113080214938675295403185329139992975166793772142271258151304442513909381604146632990073523325955307898794562533668234281782537834837658309536220932535310542771290251824724954275657240999072306410018238718367949940216067381097518632850688813066669357111649740684488924064842966278395972504667003783209650209404664470473516626530812330405192303982592057596124219993836182749544178839159646186777290994412857338708438051561329739731374132824 - 2.13935779594163842233355642709119528249521925823403989240099047419182844594163372694154617394239007198294095347993243527595277010822795133074387676608036746611809657810566130554997964611843067325866549759131456247878514062226540169569226242329355713176791248500861356587694829916732142782092779675878417336917923909520067087992398504000426262545825291783555285482826554036037196417194281024790664045908394233089603430439291179577887349850954176225940687576575164974314437660565471134851615019004153512755279693E-525*I -5.39081511188198310179590041653676667757150727264901776126586031049692562073779132715429917919427486614565251640918644362554342731912273511034790337480920545157575959751419407189957485881553540287409919517949198358871572459569827390042106968164942134433471092167805247092866858446950012358560169453254324242951479205172889272087451533604179717310455722567599516000157078206047638391319142244879511851881128877697801484429497004040128636708210057619271237107355019323892075881554403597788751760918275755820239917E-525 + 3.13779361129705433996352235915153753224586569498551077332144871703088511183161455283373809388021277936185673146410007483163521469955263521740004671025672697978420263694618384282447197413561707692134375105515702822019267627183685823159470941274428968397919893821814575180206797196864445386050599624166174601888962584604279360189098107742675960150215214839533012121012757700947588052050299597507812125553356803083527260871975234642087880690635518391842820038274727760408420183353702417271860103656763109710420701*I 5.39081511188198310179590041653676667757150727264901776126586031049692562073779132715429917919427486614565251640918644362554342731912273511034790337480920545157575959751419407189957485881553540287409919517949198358871572459569827390042106968164942134433471092167805247092866858446950012358560169453254324242951479205172889272087451533604179717310455722567599516000157078206047638391319142244879511851881128877697801484429497004040128636708210057619271237107355019323892075881554403597788751760918275755820239917E-525 - 3.13779361129705433996352235915153753224586569498551077332144871703088511183161455283373809388021277936185673146410007483163521469955263521740004671025672697978420263694618384282447197413561707692134375105515702822019267627183685823159470941274428968397919893821814575180206797196864445386050599624166174601888962584604279360189098107742675960150215214839533012121012757700947588052050299597507812125553356803083527260871975234642087880690635518391842820038274727760408420183353702417271860103656763109710420701*I]
*/


