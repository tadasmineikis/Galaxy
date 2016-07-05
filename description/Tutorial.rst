Programinio paketo diegimo procedūra
====================================

IN PROGRESS

Modelio įvedimo failų paruošimas
================================

Prieš paleidžiant modelį reikia paruošti įvesties failus:
1. Disko sukimosi kreivė. 
2. Dujų akrecijos galaktikos diske radialinis profilis.
3. Dujų akrecijos galaktikoje laiko evoliucija. 
4. SSP evoliucijos failai.

1. Disko sukimosi kreivė
~~~~~~~~~~~~~~~~~~~~~~~~

Failas prasideda neribotu komentarų kiekiu, kurie kiekvienoje naujoje
eilutėje prasideda '#' simboliu.

Failas turi turėti 2 stulpelius: 
- pirmąjame nurodyti modelio žiedų radialiniai atstumai; 
- antrąjame nurodytas modelio žiedo sukimosi greitis km/s. Skaitoma į float kintamųjų masyvą.

Pastabos:
'''''''''

-  Pirmojo sukimosi kreivės stulpelio, kuriame nurodyti atstumai modelis
   nepanaudoja, jis reikalingas tik modelio naudotojo patogumui. Modelio
   žiedų radialiniai atstumai užduodami modelio parametrų faile.

-  Failo pavyzdys:

    #r[kpc] vr[km/s]
          0.00 0.00
          0.05 1.00
          0.10 2.00
          0.15 3.00
          0.20 4.00


2. Dujų akrecijos galaktikoje laiko evoliucija
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Failas prasideda neribotu komentarų kiekiu, kurie kiekvienoje naujoje
eilutėje prasideda '#' simboliu.

Failas turi turėti 2 stulpelius: 
- pirmąjame nurodyti modelio amžiai; 
- antrąjame nurodyta į galaktiką įkritusi dujų masė Saulės masėmis. Skaitoma į float kintamųjų masyvą.

Pastabos:
'''''''''

-  Pirmasis failo stulpelis su laiko žingsnių amžiais programoje
   nenaudojamas, jis reikalingas tik modelio naudotojo patogumui.
   Modelio laiko žingsniai užduodami modelio parametrų faile.
-  Jeigu failo eilučių skaičius bus mažesnis nei modelio parametrų faile
   nurodytas laiko žingsnių kiekis, programa išluš
-  Susumuota visa į galaktiką įkrentančių dujų masė išvedama į terminalą
   paleidžiant programą

-  Failo pavyzdys

    #time[Gyr] acc[Msol]
        0.000 36943.350
        0.005 36943.350
        0.010 36943.350
        0.015 36924.575
        0.020 36905.800


3. Dujų akrecijos galaktikos diske radialinis profilis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Failas prasideda neribotu komentarų kiekiu, kurie kiekvienoje naujoje
eilutėje prasideda '#' simboliu

Failas turi turėti 2 stulpelius: 
- pirmąjame nurodyti modelio žiedų radialiniai atstumai; 
- antrąjame nurodytas radialinis dujų akrecijos profilis. Skaitoma į float kintamųjų masyvą.

Šis profilis nurodo, kokiomis proporcijomis dujos įkritusios į galaktiką bus paskirstytos
diske. Pvz.: jei modelis būtų sudarytas iš 3 žiedų su proflio vertėmis
1, 0.5, 0.1, į centrinę ląstelę įkritusių dujų būtų 2 kart daugiau nei į
2-ąją ir 10 kart daugiau nei į 3-ąją.

Pastabos.
'''''''''

-  Pirmasis failo stulpelis su radialiniais atstumais modelyje
   nenaudojamas. Jis skirtas tik programos naudotojo patogumui
-  Radialinio profilio vienetai nėra svarbūs. Programai nuskaičius failą
   profilis yra pernormuojamas taip, kad visų modelio ląstelių suma būtų
   lygi 1. Pateiktame faile naudojami vienetai Msol/pc^2, nes tuomet
   labai patogu tokį profilį interpretuoti, kaip dujų ir žvaigždžių
   paviršinių tankių profilių sumą.

Failo pavyzdys
``````````````
    #r[kpc] Tmas[Msol/pc^2]
                0.00 17.33
                0.05 15.33
                0.10 13.64
                0.15 12.45
                0.20 11.57


4. SSP (simple stellar population) evoliucijos failai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modelio veikimui reikalingi: - [gas_matrix.pegase] SSP ląstelės dujoms
grąžinama dujų masė (M_gas) ties kiekvienu laiko žingsniu ir
metalingumu. Stulpeliuose nurodyti skirtingi metalingumai, o eilutėse
skirtingi laiko žingsniai.

-  [metals_matrix.pegase] SSP ląstelės dujoms grąžinama metalų (M_gas
   x Z_gas) masė. Stulpeliuose nurodyti skirtingi metalingumai, o
   eilutėse skirtingi laiko žingsniai.

-  [ages] amžių failas, kur nurodomi SSP failo eilučių laiko žingsniai.
   Eilučių skaičius turi atitikti eilučių skaičių mgas_matrix.pegase ir
   metals_matrix.pegase failuose.

-  [zlist] metalingumų failas, kur nurodamas SSP failo stulpelių
   metalingumai. Eilučių skaičius turi atitikti stulpelių skaičių
   mgas_matrix.pegase ir metals_matrix.pegase failuose.

-  [phot_u.flux] fotometrijos failas. Jis nėra būtinas modelio
   veikimui, tačiau reikalingas norint išvesti fotometrijos rezultatus.
   Stulpeliuose nurodomas sunormuotas į 1 Msol fliuksas ties skirtingais
   metalingumais, o eilutėse ties skirtingais laiko žingsniais.

Šių failų pavyzdžiai pateikiami **description/** direktorijoje. Šių
failų [phot-u.flux, metals-matrix.pegase, gas-matrix.pegase] paruošimo
skriptai pateikiami **tools/** direktorijoje

Modelio parametrų failų paruošimas
==================================

1. Sintetinės fotometrijos parametrų failas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Programa skaitydama parametų failą eilutes prasidedančias simboliu #
praleidžia kaip komentarus. Nesant šio simbolio eilutė interpretuojama
pagal nurodytą parametro vardą.

1. ISO num_of_files :

   -  num_of_files integer tipo kintamas, po šio parametro programa
      sekančioje eilutėje ieškos nurodyto kelio iki kiekvieno izochronos
      failo.

2. IMF sampling_method integration num_of_rows num_of_cols :

   -  smapling_method : string tipo parametras leidžiantis pasirinkti
      skirtingus pradinių masių funkcijos generavimo metodus,
      dažniausiai naudojama opcija - 'stochastic'
   -  integration : int tipo parametras, kuris buvo naudingas
      generuojant spiečių integruotą fotometriją, leisdavo pasirinkti,
      kurią dalį žvaigždžių šviesio sumuoti stochastiškai, o kurią
      analitiškai. Generuojant CMD reikia nustayti vertę '0'
   -  num_of_rows : int tipo kintamasis, nurodantis kiek intervalų
      turės pradinių masių funkcija.
   -  num_of_cols : int tipo kintamasis, nurodantis kiek stulpelių
      turės pradinės masių funkcijos aprašymas. Pirmąjame stulpelyje
      nurodomas masių intervalo apatinė masė, antrąjame - viršutinė,
      trečiąjame - laipsninis rodiklis, nurodomas `alpha` forma.

3. SEED seed:

   -  seed : int tipo atsitiktinių skaičių generatoriaus seed'as. Ši
      eilutė užpildoma automatiškai leidžiant python script'ą

4. BINNARY_FRACTION binnary_fraction

   -  binnary_fraction : float tipo parametras nurodantis dvinarių
      sistemų dalį žvaigždžių populiacijoje

5. LIMIT filter_num magnitude_cut

   -  filter_num : int tipo parametras atitinkantis parametrą, pagal
      kurį apriboti išvedamų žvaigždžių kiekį. filter_num atitinka
      stulpelį, iš užduoto izochronos failo.
   -  magnitude_cut : float tipo kintamasis pagal kurį apribojamas
      išvedimas.

6. THREADS num_of_threads

   -  num_of_threads : int tipo kitamasis, generuojant CMD žvaigždes
      veikia tik su opcija '1'

7. GALEMO_RESULTS

   -  nurodo nuskaitomo failo tipą, šis eilutė automatiškai užpildoma
      leidžiant python script'ą

8. OUT

-  nurodo išvedimo failo vardą, šis eilutė automatiškai užpildoma
   leidžiant python script'ą

Failo pavyzdys
``````````````

    #============================================
    iso 3
    iso_bank/UBV/parsec1.2s_dense/icz0001.dat
    iso_bank/UBV/parsec1.2s_dense/icz0002.dat
    iso_bank/UBV/parsec1.2s_dense/icz0003.dat
    #********************************************
    imf stochastic 0 4 3
    0.01 0.08 -0.3
    0.08 0.5 -1.3
    0.5 1.0 -2.3
    1.0 120 -2.7
    #********************************************
    SEED
    BINARY_FRACTION 0.
    #============================================
    limit 4 4
    #********************************************
    threads 1
    GALEMO_RESULTS
    #
    OUT
    #********************************************


2. Modelio parametrų failas
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1.  TIMESTEP timestep:

    -  timestep [Myr]: int tipo kintamasis nurodantis modelio
       integravimo laiko žingsnį

2.  GALAXY_AGE galaxy_age:

    -  galaxy_age [Myr]: int tipo kintamasis nurodantis modelio
       integravimo laiką

3.  SN_TIMESCALE sn_timescale:

    -  sn_timescale [Myr]: int tipo kitamasis nurodantis laiko
       intervalą, kurį ląstelėje bus palaikomas mažas dujų tankis,
       praturtintų metalais dujų praradimas iš galaktikos ir visiškai
       negalima žvaigždėdara

4.  SEED seed:

    -  seed : int tipo kitamasis nurodantis atsitiktinių skaičių
       generatoriaus pradinį 'seed'. Naudojamas norint atkartoti
       modeliavimo rezultatus

5.  DISTANCE distance:

    -  distance [Mpc]: float tipo kintamasis nurodantis modeliuojamos
       galaktikos atstumą. Naudojams mag/arcsec^2 paviršinio šviesio
       išvedimui ouput failuose

6.  GRID_SIZE grid_size grid_buffer

    -  grid_size: int tipo kintamasis nurodantis disko žiedų skaičių
    -  grid_buffer: int tipo kintamasis nurodantis kiek žiedų pridėti
       prie 'grid_size' kraštiniams efektams modelyje panaikinti,
       priklausomai nuo grid_size galima pridėti nuo0 iki 20%
       grid_size vertės

7.  CELL_SIZE cell_size

    -  cell_size [pc]: int tipo kintamasis nurodantis ląstelės fizinius
       matmenis

8.  GAS_DIFFUSION gas_diff

    -  gas_diff [Myr]: float tipo kintamasis nurodantis dujų tankio
       netolygumų tarp ląstelių išslyginimo laiko skalę

9.  STELLAR_DIFFUSION stell_diff mass_treshold

    -  stellar diff: float tipo kintamasis nurodantis kuri žvaigždžių
       masės dalis per laiko žingsnį prarandama/gaunama į/iš kaimyninių
       ląstelių
    -  mass_treshold [Msol]: float tipo kintamasis nurodantis minimalią
       apsikeitimo žvaigždėmis masę tarp *visų kaimyninių ląstelių*, jei
       suminė apsikeitimo masė yra mažesnė už šį parametrą, apsikeitimas
       žvaigždėmis tarp ląstelių nevyksta.

10. SFE epsilon

    -  epsilon : float arba string tipo kintamas. Jei tai float kintamasis,
       jis nurodo žvaigždėdaros efektyvumą įvykus žybsniui ląstelėje
       (pvz.:0.05 atitinka 5% efektyvumą); jei tai string kintamasis,
       jis interpretuojamas kaip failo vardas, kuriame nurodytos sfe
       vertės kiekvienam modelio laiko žingsniui (failas vieno stulpelio).

11. TRIGGERED trig trig_time

    -  trig: float arba string tipo kintamasis. Jei tai float
       kintamasis, jis nurodo indukuotos žvaigždėdaros tikimybę
       kaimyninėse ląstelėse (kiekvienai ląstelei atskirai, suminė
       tikimybė būtų kaimynių skaičius x trig). Jei tai string
       kintamasis, jis nurodo failą kur trig vertės nurodytos kiekvienam
       modelio laiko žingsniui (failas vieno stulpelio).
    -  trig_time [time steps] : int tipo kintamasis kuris nurodo kuriuo
       laiko žingsniu po žvaigždėdaros žybsnio vyksta kaimyninių
       ląstelių indukavimas (pvz.:1 - sekančiu laiko žingsniu, 2 - po
       dviejų laiko žingsnių)

12. SFE_POW alpha

    -  alpha: float arba string kintamasis. Jei tai float kintamasis,
       jis nurodo žvaigždėdaros efektyvumo priklausomybės nuo dujų
       tankio laispninį rodiklį (pvz.: 0 -- nepriklauso, 1 -- tiesiškai,
       2 -- kvadratinė priklausomybė); jei tai string kintamasis, jis
       nurodo failą, kur šis parametras nurodytas kiekvienam modelio
       laiko žingsniui (failas vieno stulpelio).

13. NRM_SFE nrm_sfe

    -  nrm_sfe [Msol/pc^2] - float tipo kintamasis, kuris nurodo ties
       kokiu dujų tankiu normuojamas žvaigždėdaros efektyvumo dėsnis.
       Pvz.: jei nrm_sfe=10 Msol/pc^2, tai reiškia, kad ląstelės dujų
       pavišiniui tankiui esant 10 Msol/pc^2, įvykus žvaigždėdaros
       žybsniui SFE bus lygus epsilon

14. MINIMUM_SFE minimum_sfe

    -  minimum_sfe : float tipo kintamsis nurodantis minimalų
       žvaigždėdaros efektyvumą žvaigždėdaros žybsnio metu, jei
       žvaigždėdaros efektyvumas dėl dujų tankio pagal užduotą dėsnį
       būtų mažesnis.

15. MAXIMUM_SFE maximum_sfe

    -  maximum_sfe : float tipo kintamasis nurodantis maksimalų
       žvaigždėdaros efektyvumą žvaigždėdaros žybsnio metu, jei pagal
       užduotą žvaigždėdaros efektyvumas dėl dujų tankio pagal užduotą
       dėsnį būtų didesnis

16. TRIGG_MASS minimum_mass trigg_mass

    -  minimum_mass [Msol] : float tipo kintamasis, nurodantis, kokios
       mažiausios susiformavusios žvaigždžių masės žvaigždėdaros žybsnis
       gali vykti ląstelėje. Jei žvaigždžių masė gimusi ląstelėje būtų
       mažesnė už šį dydi, žvaigždėdaros žybsnis nevyksta
    -  trigg_mass [Msol] : float tipo kintamasis, nurodantis kokios
       masės turi būti susiformavusi žvaigždžių populiacija po
       žvaigždėdaros žybsnio, kag galėtų indukuoti kaimynines ląsteles

17. GAS2_SPONT gas2_spont

    -  gas2_spont : float tipo kintamsis nurodantis kiek spotaninės
       žvaigždėdaros žybsnių vidutiniškai vyks per vieną laiko žingsnį
       modelyje. Tikimybė vykti žvaigždėdarai ląstelėje yra proporcinga
       dujų paviršinio tankio kvadratui.

18. GAS_SFR_TRESHOLD gas_sfr_treshold

    -  gas_sfr_treshold [Msol/pc^2] : float tipo kintamasis nurodantis
       ribinį dujų paviršinį tanki žvaigždėdarai vykti. Jei ląstelės
       paviršinis dujų tankis yra mžesnis už ribinį, žvaigždėdaros
       tikimybė lastelėje tiesiškai mažinama.

19. OUTFLOW eta critical_velocity

    -  eta: float tipo parametras nurodantis kuri dalis metalais
       praturtintų gražinamų visų evoliucionuojančių žvaigždžių
       populiacijų pararandama.
    -  critical_velocity : float tipo parametras nurodantis kokios
       energijos turi būti žvaigždėdaros žybsnis, kad išmestų dujas į
       tarp galaktinę erdvę. Parametro fizikinė prasmė - besiplečinčio
       superburbulo greitis prieš 'greitėjimą'. 1 atitinka 20 km/s,
       pagal straipsnio autorius, pakankamas greitis išmesti dujas iš
       Paukščių Tako disko.

20. ACC_METALLICITY acc_metallicity

    -  acc_emtalicity : float tipo parametras nurodantis akrecijos dujų
       metalingumą diske užduodantį failą.

21. ACRETION acc_rad_file acc_time_file

    -  acc_rad_file : string tipo parametras nurodantis akrecijos
       radialinę priklausomybę diske užduodantį failą.
    -  acc_time_file : string tipo kintamasis nurodantis akrecijos
       laiko priklausomybės galaktikoje failą.

22. OUTPUT num_of_types types[] num_of_times times[]

    -  num_of_types : int tipo kintamasis kuris nurodo kiek bus
       išvedama failų tipų
    -  types[] : string tipo kintamųjų grupė, kuri nurodo išvedimo
       tipus, galimos opcijos: 0d, 1d, 2d
    -  num_of_times : int tipo kintamasis nurodantis kiek bus išvedimo
       laikų
    -  times[] [Myr]: int tipo kintamųjų grupė, nurodanti kuriais
       modelio amžiais išvesti nurodytus failus

23. NUM_OF_THREADS num_of_threads

    -  num_of_threads : int tipo kintamasis nurodantis modeliui
       skaičiuoti naudojamų gijų skaičių.

24. CMD_LIMIT mass_limit age_limit

    -  mass_limit [Msol] : float tipo kintamasis nurodantis mažiausią
       išvedamos SSP populiacijos masę
    -  age_limit [Myr] : float tipo kintamasis nurodantis išvedamos
       populiacijos maksimlų amžių.

25. ROTATION rotation_curve_file

    -  rotation_curve_file : string tipo kintamasis nurodantis disko
       suskimosi kreivę užduodančio failą vardą

26. PEGASE mgas_file zgas_file

    -  mgas_file : string tipo kintamasis su failo vardu, kuriame yra
       SSP prarandamų dujų kiekiai kiekvienam SSP metalingumui ir amžiui
    -  zgas_file : string tipo kintamsis su failo vardu, kuriame yra
       SSP prarandamų metalų (Z_gas x M_gas) kiekiai kiekvienam SSP
       amžiui ir metalingumui.

27. SSP ages zlist

    -  ages : string tipo kintamasis nurodantis failo vardą su modelio
       laiko žingsnių amžiais
    -  zlist : string tipo kintamasis nurodantis failą su SSP
       metalingumo žingsniais

28. PHOTOMETRY num_of_filters

    -  num_of_filters : int tipo kintamsis nurodantis kiek SSP filtrų
       bus naudojama išvedime
    -  sekančiose eilutėse nurodoma filter_file filter_name porų kiek
       užduota parametru num_of_filters.

Modelio paleidimas
==================

paleidimas naudojant python pagalbinį skriptą
---------------------------------------------

Pasitikirnam, ar visi reikalingi failai yra direktorijoje, t.y.: 
- akreijos radialinio profilio failas 
- akrecijos laiko evoliucijos failas
- sukimosi kreivė 
- modelio parametrų failas 
- žvaigždžių sintetinės fotometrijos parametrų failas, pavadintas "template"

RUN direktorijoje paleidžiam python script'ą su modelio parametrų failu ir iteracijų kiekiu kaip argumentu: 

./galemo.py modelio_parametrų_failas iteracijų_skaičius

Modelio išvesties duomenų failai
================================

Modelio išvedimo failai yra sugrupuoti į 4 kategorijas: 
- 0d failuose išvedamos globalūs visos glaktikos fizikiniai parametrai; 
- 1d failuose išvedami suvidurkinti azimutiškai vienos ląstelės pločio žieduose galaktikos parametrai; 
- 2d failuose išvedami visų disko ląstelių fizikiniai parametrai;
- cmd failuose išvedami sintetiniai žvaigždžių fotometrijos katalogai.

0d išvedimai failo stulpelių paaiškinimai:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  t [Myr]: laiko žingsnis
-  SP_E: spontaninės žvaigždėdaros žybsniai galaktikoje per laiko
   žingsnį
-  TR_E: indukuotos žvaigždėdaros žybsniai galaktikoje per laiko
   žingsnį
-  ACC [Msol/yr]: dujų akrecijos sparta galaktikoje
-  ST_GAS_ACC [Msol/yr]: evoliucionuojančių SSP gražinamų dujų sparta
-  OTFL [Msol/yr]: prarandamų praturtintų sunkiaisiais metalais dujų
   sparta
-  STARS [Msol]: žvaigždžių masė modelio diske
-  GAS [Msol]: dujų masė modelio diske
-  ZGAS: dujų metalingumas modelio diske
-  TSFR [Msol/yr]: žvaigždėdaros sparta galaktikoje
-  ModelAge [Myr]: modelio galautinis amžius
-  iteration: modelio iteracijos indeksas

1d išvedimo failo stulpelių paaiškinimai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  r [kpc] : žiedo vidurio radialinis atstumas iki disko centro
-  mgas [Msol/pc^2] : disko žiedo dujų vidutinis paviršinis tankis
-  zgas : disko žiedo dujų medianinis metalingumas
-  mstr [Msol/pc^2] : disko žiedo žvaigždžių vidutinis paviršinis tankis
-  SF_events : suminis žvaigždėdaros žybsnių skaičius disko žiede per visą modeliavimo laiką
-  SP_events : suminis spontaninės žvaigždėdaros žybsnių skaičius žiede per visą modeliavimo laiką
-  SFR [Msol/pc^2/Myr] : žvaigždėdaros spartos paviršinis tankis diske žiede
-  SFR100 [Msol/pc^2/ (10 x timesteps) Myr] : žvaigždėdaros spartos per paskutinius 10 laiko žingnsių paviršinis tankis
-  ACTIVE : jeigu aktyvuota, skaičiuoja suminį ląstelių galinčių formuoti žvaigždes skaičių
-  Tgas [Msol/pc^2] : iki išvedimo laiko į disko ląsteles įkritusių dujų 
   pavišinis tankis įkritusių (jei nevykst ažvaigždėdara, turi būti
   lygus mgas stulpeliui, arba apytiksliai lygus mgas+mstr)
-  Ogas-tot [Msol] : vidutinė iš 1-os ląstelės prarasta (outflow)
   praturtinų metalais dujų masė per modeliavimo laiką
-  Ogas_cur [Msol/pc^2/ (10 x timesteps) Myr] : vidutinė iš 1-os
   ląstelės prarandama praturtintų dujų masė per 10 laiko žingsnių
-  Ometals-tot [Msol] : vidutinė iš 1-os ląstelės prarasta metalų masė
   (Mgas\*Zgas) per modeliavimo laiką
-  Ometals-cur [Msol/pc^2/(10 x timesteps) Myr] : vidutinė iš vienos
   ląstelės prarandama metalų masė per 10 laiko žingnsių
-  ModelAge [Myr]: išvedimo amžius
-  iteration : modelio iteracijos indeksas

2d išvedimo failo stulpelių paaiškinimai:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  r [kpc] : ląstelės vidurio radialinis atstumas iki disko centro
-  x [kpc] : ląstelės centro koordinatė xy plokštumoje
-  y [kpc] : ląstelės centro koordintaė xy plokštumoje
-  ref_t [timesteps] : laiko žingsniai praėjo po paskutinio žvaigždėdaros žybsnio ląstelėje
-  sfr_t [ratio] : dujų masės santykis ląstelėje lyginant su kritiniu tankiu žvaigždėdarai
-  sp_buff : spotaninės žvaigždėdaros parametras naudotas 'debuggininimui'
-  mgas [Msol/pc^2] : paviršinis dijų tankis ląstelėjė
-  zgas : dujų metalingumas ląstelėje
-  mstr [Msol/pc^2]: žvaigždžių paviršinis tankis ląstelėje
-  sfe [] : pakutinio žvaigždėdaros žybsnio efektyvumas (M_str/(M_gas+M_str))
-  last_mstr [Msol] : paskutinio žvaigždėdaros žybsnio metu susiformavusios žvaigždžių populiacijos masė
-  TVEL : dujų praradimo parametras, rodantis besiplečiančio burbulo greičio santyki su kritiniu
-  rho : dujų dalelių tankis į cm^3 disko plokštumoje, darant prielaidą, jog dujos diske Z kryptimi pasiskirsčiusios eksponentiškai
-  N0 : outflow paramertas iš Baumgartner V., Breitschwerdt, D. 2013, AA, 557, A140
-  ModelAge [Myr]: išvedimo amžius
-  iteration : modelio iteracijos indeksas

cmd išvedimo failo stulpelių paaiškinimai:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  r [kpc] : ląstelės centro, kurioje yra žvaigždė, radialinis atstumas nuo galaktikos centro
-  a [rad] : ląstelės centro, kurioje yra žvaigždės, azimutinė koordiantė išreikštos radianais
-  x [kpc] : koordinatė xy plokštumoje, su triukšmu paskirstančiu žvaigždes tolygiai ląstelės ribose
-  y [kpc] : koordinatė xy plokštumoje, su triukšmu paskirstančiu žvaigždes tolygiai ląstelės ribose
-  age [log_10(t/yr)] : žvaigždei priskirtas artimiausios izochronos amžius
-  a_age [log_10(t/yr)] : tikrasis žvaigždės amžius
-  z : žvaigždei priskirtas artimiausisas izochronos metalingumas
-  a_z : tikrasis žvaigždės metalingumas
-  mass_i_p [log_10(m/Msol)] : žvaigždės / pirminės dvianrės
   žvaigždės masė
-  mass_i_s [log_10(m/Msol)] :  antrinės dvinarės žvaigždės masė
-  p_X [mag] :  pirminės dvinarės sistemos žvaigždės
   fotometrinis ryškis
-  s_X [mag] :  antrinės dvinarės sistemos žvaigždės fotometrinis ryškis
-  o_X [mag] : suminis dvinarės sistemos fotometrinis ryškis
-  ssp_index : SSP indeksas, kad būtų galima identifikuoti žvaigždes iš tos pačios SSP
-  ModelAge [Myr] : modelio amžius išvedimo metu
-  iteration : modelio iteracijos indeksas

Modelio išvesties grafiniai failai
==================================

Skaičiuojant modelius python script'o 'galemo' pagalba, baigus skaičiavimus
sukuriame direktorija Dir_+'parametrų failo vardas', kurioje sukuriami
pradinės analizės diagnostiniai paveiksliukai.

0d tipo paveiksliukai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sukuriama po vieną kiekvienai iteracijai. Paveiksliuke pavaizduoti
6 grafikai (3 eilutės, 2 stulpeliai), visų grafikų Y ašis logaritminė:

1-1 (nuo viršaus kairės) - Žvaigždėdaros istorija galaktikoje
nuo laiko. Kiekvienas taškas grafike atitinka žvaigždėdaros spartą kiekvienu
laiko žingsniu.

1-2 - dujų akrecijos sparta galaktikoje.

2-1 - spontaninių žvaigždėdaros žybsnių evoliucija galaktikoje. Kiekvienas
taškas atitinka spontatninės žvaigždėdaros įvykių skaičių per laiko žingsnį.

2-2 - indukukuotos žvaigždėdaros žybsnių evoliucija galaktikoje. Kiekvienas
taškas atitinka indukuotos žvaigždėdaros žybsnių skaičių per laiko žingsnį.

3-1 - sunkiaisiai elementais praturtintų dujų praradimo evoliucija.

3-2 - dujų, grąžinamų evoliucionuojančių SSP populiacijų į galaktikos tarpžvaigždinę terpę
sparta

1d tipo paveiksliukai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sukuriama po vieną kiekvienam išvedimo laikui ir iteracijai. Paveiksliuke pavaizduota
12 grafikų (4 eilutės, 3 stulpeliai), visų grafikų Y ašis logaritminė:

1-1 (nuo viršaus kairės) - dujų paviršinio tankio profilis glaktikoje.

1-2 - žvaigždžių paviršinio tankio profilis galaktikoje.

1-3 - dujų metalingumo profilis galaktikoje (Y ašyje pavaizduota log(Z)).

2-1 - žvaigždėdaros spartos (per laiko žingsnį) profilis galaktikoje.

2-2 - žvaigždėdaros spartos (per 10 laiko žingsnių) profilis galaktikoje.

2-3 - vidutinė iš ląstelės žiede dėl outflow prarandama suminė dujų masė.

3-1 - vidutinė iš ląstelės žiede dėl outflow prarandama dujų masė per 10 laiko žingsnių.

3-2 - vidutinė iš ląstelės žiede dėl outflow prarandama suminė metalų masė.

3-3 - vidutinė iš ląstelės žiede dėl outflow prarandama metalų masė per 10 laiko žingsnių.

4-1 - suminis žvaigždėdaros įvykių (spontaninių ir indukuotų) profilis galaktikoje

4-2 - suminis spontaninių žvaigždėdaros įvykių profilis galaktikoje

4-3 - Dėl akrecijos įkrituios suminis dujų masės tankio profilis galaktikoje. Realus masės tankio profilis
sutaptų su šiuo tik išjungus žvaigždėdarą galaktikoje. Dėl žvaigždžių dispersijos,
dujų judėjimų tarp ląstelių ir dujų outflow realus profilis atrodo kitaip.


2d tipo paveiksliukai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sukuriama po vieną kiekvienam išvedimo laikui ir iteracijai. Paveiksliuke pavaizduoti
8 grafikai (4 eilutės, 2 stulpeliai), visų grafikų Y ašis logaritminė:

1-1 (nuo viršaus kairės) - dujų paviršinio tankio profilis glaktikoje.

1-2 - žvaigždžių paviršinio tankio profilis galaktikoje.

2-1 - dujų metalingumo profilis galaktikoje (Y ašyje pavaizduota log(Z)).

2-2 - paskutinio žvaigždėdaros žybsnio metu buvusio žvaigždėdaros efektyvumo
(Mstr/(Mgas+Mstr)) profilis galaktikoje.

3-1 - paskutinio žvaigždėdaros žybsnio metu susiformavusiios žvaigždžių masės profilis galaktikoje.

3-2 - paskutinio žvaigždėdaros žybsnio ląstelėje amžiaus (laiko žingsniais) profilis galaktikoje.

4-1 - dujų tankio ląstelėje santykio su kritiniu dujų tankiu profilis galaktikoje.

4-2 - diagnostinis outflow prametro TVEL profilis galaktiko:
- -9.99  - susiformavo 1 supernova arba mažiau
- -99.99 - outflow išjungtas
-  9999  - susiformavo daugiau nei 500 supernovų
-  2     - supernovų kiekis pakankamas, kad vyktų outflow
- -1     - supernovų kiekis nepakankamas, kad vyktų outflow




cmd tipo paveiksliukai
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sukuriama po vieną kiekvienam išvedimo laikui ir iteracijai.
Paveiksliuką sudaro 3 grafikai:

- 1-1 (nuo viršaus kairės) CMD I-BI filtruose
 (jei šių filtrų išvedimo failuose nėra, reikės paredaguoti script'ą)
- 2-1 - žvaigždėdaros spartos evoliucija. Juodi taškai žymi žvaigždėdaros spartą
kiekvienam laiko žingsniui, histograma (juoda linija) vidurkį, raudona linija - dujų akrecijos spartą

- 2-2 - juoda linija žymi dujų metalingumo evoliuciją, 
raudona - žvaigždžių masės evoliuciją

