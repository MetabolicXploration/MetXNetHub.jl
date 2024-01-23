## --------------------------------------------------------
# Utils
function _extract_rxnid(line::AbstractString)
    rxn_id_reg = r"v[A-Za-z0-9]+\s*\(.+\)"
    m = match(rxn_id_reg, line)
    return isnothing(m) ? "" : m.match
end

function _extract_metid(word::AbstractString)
    met_id_reg = r"[A-Za-z][A-Za-z0-9]+"
    m = match(met_id_reg, word)
    return isnothing(m) ? "" : m.match
end

function _extract_stoi_coe(word::AbstractString)
    reg = r"^[0-9]+$"
    m = match(reg, word)
    return isnothing(m) ? "" : m.match
end

function _extract_arrow(word::AbstractString)
    reg = r"^→$"
    m = match(reg, word)
    return isnothing(m) ? "" : m.match
end


function _is_match(reg::Regex, w::String)
    m = match(reg, w)
    return !isnothing(m)
end

# ------------------------------------------------------------------
# The string was extracted from 
# Massucci et al., “Energy Metabolism and Glutamate-Glutamine Cycle in the Brain.”
# Supplementary materials pdf

const _MASSUCCI_NETWORK_RAW_STR = """
1 vGLC(c) → GLCc
2 vGLC(c → a) GLCc → GLCa
3 vGLC(c → e) GLCc → GLCe
4 vGLC(e → a) GLCe → GLCa
5 vGLC(e → n) GLCe → GLCn
6 vO2(c) → O2c
7 vO2(c → a) O2c → O2a
8 vO2(c → n) O2c → O2n
9 vLAC(c) LACc →
10 vLAC(e → c) LACe → LACc
11 vLAC(a → c) LACa → LACc
12 vLAC(e → a) LACe → LACa
13 vLAC(a → e) LACa → LACe
14 vLAC(e → n) LACe → LACn
15 vLAC(n → e) LACn → LACe
16 vCKf (n) ATPn + Crn → ADPn + PCrn
17 vCKr(n) ADPn + PCrn → ATPn + Crn
18 vCKf (a) ATPa + Cra → ADPa + PCra
19 vCKr(a) ADPa + PCra → ATPa + Cra
20 vAKf (n) 2 ADPn → AMPn + ATPn
21 vAKr(n) AMPn + ATPn → 2 ADPn
22 vAKf (a) 2 ADPa → AMPa + ATPa
23 vAKr(a) AMPa + ATPa → 2 ADPa
24 vPPP(n) 3 G6Pn + 6 NADPn → 2 F6Pn + GAPn + 6 NADPHn
25 vPPP(a) 3 G6Pa + 6 NADPa → 2 F6Pa + GAPa + 6 NADPHa
26 vNAK(n) ATPn + 2 Ke + 3 Nan → ADPn + 2 Kn + 3 Nae
27 vNAK(a) ATPa + 2 Ke + 3 Naa → ADPa + 2 Ka + 3 Nae
28 vNKCC(a) Ke + Nae → Ka + Naa
29 vNa(n) Nae → Nan
30 vK(n) Kn → Ke
31 vGR1(n) GSSGn + NADPHn → GSHn + NADPn
32 vGR2(n) GSSGn + NADPHnm → GSHn + NADPnm
33 vDHAR(n) DHAn + GSHn → ASCn + GSSGn
34 vAPX(n) ASCn + ROSn → DHAn
35 vGR1(a) GSSGa + NADPHa → GSHa + NADPa
36 vGR2(a) GSSGa + NADPHam → GSHa + NADPam
37 vDHAR(a) DHAa + GSHa → ASCa + GSSGa
38 vAPX(a) ASCa + ROSa → DHAa
39 vGSH(a) → GSHa
40 vGSH(a → e) GSHa → GSHe
41 vGSH(e → n) GSHe → GSHn
42 vDHA(n → e) DHAn → DHAe
43 vDHA(e → n) DHAe → DHAn
44 vDHA(e → a) DHAe → DHAa
45 vDHA(a → e) DHAa → DHAe
46 vASC(a → e) ASCa → ASCe
47 vASC(e → a) ASCe → ASCa
48 vASC(e → n) ASCe + 2 Nae → ASCn + 2 Nan
49 vGLU(e → a) GLUe + Ka + 3 Nae → GLUac + Ke + 3 Naa
50 vGS(a) ATPa + GLUac → ADPa + GLNa
51 vGLN(a → n) GLNa → GLNn
52 vPAG(n) GLNn → GLUnc
53 vGLU(n) ATPn + GLUnc → ADPn + GLUnv
54 vNT(n → e) GLUnv → GLUe
55 vHK(n) ATPn + GLCn → ADPn + G6Pn
56 vPFK(n) ATPn + G6Pn → ADPn + 2 GAPn
57 vGAPDHf (n) GAPn + NADnc → BPGn + NADHnc
58 vGAPDHr(n) BPGn + NADHnc → GAPn + NADnc
59 vPGKf (n) ADPn + BPGn → ATPn + PEPn
60 vPGKr(n) ATPn + PEPn → ADPn + BPGn
61 vPK(n) ADPn + PEPn → ATPn + PYRnc
62 vLDHf (n) NADHnc + PYRnc → LACn + NADnc
63 vLDHr(n) LACn + NADnc → NADHnc + PYRnc
64 vHK(a) ATPa + GLCa → ADPa + G6Pa
65 vPFK(a) ATPa + G6Pa → ADPa + 2 GAPa
66 vGAPDHf (a) GAPa + NADac → BPGa + NADHac
67 vGAPDHr(a) BPGa + NADHac → GAPa + NADac
68 vPGKf (a) ADPa + BPGa → ATPa + PEPa
69 vPGKr(a) ATPa + PEPa → ADPa + BPGa
70 vPK(a) ADPa + PEPa → ATPa + PYRac
71 vLDHf (a) NADHac + PYRac → LACa + NADac
72 vLDHr(a) LACa + NADac → NADHac + PYRac
73 vPYRf (n) PYRnc → PYRnm
74 vPYRr(n) PYRnm → PYRnc
75 vPYRf (a) PYRac → PYRam
76 vPYRr(a) PYRam → PYRac
77 vOP(n) 5 ADPn + 2 NADHnm + O2n → 5 ATPn + 2 NADnm + 0.01 ROSn
78 vOP(a) 5 ADPa + 2 NADHam + O2a → 5 ATPa + 2 NADam + 0.01 ROSa
79 vNAD(n) NADHnm + NADPnm → NADPHnm + NADnm
80 vNAD(a) NADHam + NADPam → NADPHam + NADam    
81 vPC(a) ATPa + PYRam → ADPa + OAAam
82 vcME(n) MALnc + NADPn → NADPHnc + PYRnc
83 vmME(n) MALnm + NADPnm → NADPHnm + PYRnm
84 vcME(a) MALac + NADPa → NADPHac + PYRac
85 vmME(a) MALam + NADPam → NADPHam + PYRam
86 vPDH(n) CoAn + NADnm + PYRnm → ACoAn + NADHnm
87 vCS(n) ACoAn + OAAnm → CITn + CoAn
88 vIDH1(n) CITn + NADnm → AKGnm + NADHnm
89 vIDH2(n) CITn + NADPnm → AKGnm + NADPHnm
90 vIDH3(n) CITn + NADPn → AKGnc + NADPHn
91 vAKGDH(n) AKGnm + CoAn + NADnm → NADHnm + SCoAn
92 vSCoATKf (n) ADPn + SCoAn → ATPn + SUCn
93 vSCoATKr(n) ATPn + SUCn → ADPn + SCoAn
94 vSDHf (n) 1.50 ADPn + 0.10 O2n + 3 SUCn → 1.50 ATPn + 3 FUMn
95 vSDHr(n) 1.50 ATPn + 3 FUMn → 1.50 ADPn + 0.10 O2n + 3 SUCn
96 vFUMf (n) FUMn → MALnm
97 vFUMr(n) MALnm → FUMn
98 vmMDH(n) MALnm + NADnm → NADHnm + OAAnm
99 vPDH(a) CoAa + NADam + PYRam → ACoAa + NADHam
100 vCS(a) ACoAa + OAAam → CITa + CoAa
101 vIDH1(a) CITa + NADam → AKGam + NADHam
102 vIDH2(a) CITa + NADPam → AKGam + NADPHam
103 vIDH3(a) CITa + NADPa → AKGac + NADPHa
104 vAKGDH(a) AKGam + CoAa + NADam → NADHam + SCoAa
105 vSCoATKf (a) ADPa + SCoAa → ATPa + SUCa
106 vSCoATKr(a) ATPa + SUCa → ADPa + SCoAa
107 vSDHf (a) 1.50 ADPa + 0.10 O2a + 3 SUCa → 1.50 ATPa + 3 FUMa
108 vSDHr(a) 1.50 ATPa + 3 FUMa → 1.50 ADPa + 0.10 O2a + 3 SUCa
109 vFUMf (a) FUMa → MALam
110 vFUMr(a) MALam → FUMa
111 vmMDH(a) MALam + NADam → NADHam + OAAam
112 vG3PS(n) 1.50 ADPn + 3 NADHnc + 0.10 O2n → 1.50 ATPn + 3 NADnc
113 vG3PS(a) 1.50 ADPa + 3 NADHac + 0.10 O2a → 1.50 ATPa + 3 NADac
114 vcMDH(n) NADHnc + OAAnc → MALnc + NADnc
115 vOGC(n) AKGnm + MALnc → AKGnc + MALnm
116 vAGC(n) ASPnm + GLUnc → ASPnc + GLUnm
117 vcMDH(a) NADHac + OAAac → MALac + NADac
118 vOGC(a) AKGam + MALac → AKGac + MALam
119 vAGC(a) ASPam + GLUac → ASPac + GLUam
120 vcAATf (n) AKGnc + ASPnc → GLUnc + OAAnc
121 vcAATr(n) GLUnc + OAAnc → AKGnc + ASPnc
122 vmAATf (n) AKGnm + ASPnm → GLUnm + OAAnm
123 vmAATr(n) GLUnm + OAAnm → AKGnm + ASPnm
124 vcAATf (a) AKGac + ASPac → GLUac + OAAac
125 vcAATr(a) GLUac + OAAac → AKGac + ASPac
126 vmAATf (a) AKGam + ASPam → GLUam + OAAam
127 vmAATr(a) GLUam + OAAam → AKGam + ASPam
128 vASPf (n → a) ASPnc → ASPac
129 vASPr(n → a) ASPac → ASPnc
130 vASPf (n) ASPnc → ASPnm
131 vASPr(n) ASPnm → ASPnc
132 vASPf (a) ASPac → ASPam
133 vASPr(a) ASPam → ASPac
134 vGDHf (n) GLUnm + NADPnm → AKGnm + NADPHnm
135 vGDHr(n) AKGnm + NADPHnm → GLUnm + NADPnm
136 vGDHf (a) GLUam + NADPam → AKGam + NADPHam
137 vGDHr(a) AKGam + NADPHam → GLUam + NADPam
138 vATP(n) ATPn → ADPn
139 vATP(a) ATPa → ADPa
"""
    
# ------------------------------------------------------------------
# function massucci_network(T = CoreModel)
function _Massucci2013_builder()
    
    # Parse
    dig = filter(!isempty, strip.(split(_MASSUCCI_NETWORK_RAW_STR, "\n"; keepempty = false)))

    rxn_ids = _extract_rxnid.(dig)
    
    dig = replace.(dig, rxn_ids .=> "")
    dig = replace.(dig, r"^[0-9]+" .=> "")
    dig = strip.(dig)

    # Cosmetics
    rxn_ids .= replace.(rxn_ids, r"\s" => "")
    
    S_dict = Dict()
    met_ids = []
    stoi = nothing
    sense = -1
    for (rxni, (rxn_id, met_line)) in enumerate(zip(rxn_ids, dig))
        words = strip.(split(met_line; keepempty = false))
        # parse
        for w in words
            # test stoi coe
            stoi_str = _extract_stoi_coe(w)
            if !isempty(stoi_str)
                stoi = parse(Int, stoi_str)
                continue
            end

            # test arrow
            arrow_str = _extract_arrow(w)
            if !isempty(arrow_str)
                sense *= -1
                continue
            end

            # test metid 
            met_id = _extract_metid(w)
            if !isempty(met_id)
                # I found a metabolite
                stoi = isnothing(stoi) ? 1 * sense : stoi * sense

                # store
                S_dict[(met_id, rxn_id)] = stoi
                push!(met_ids, met_id)
                
                # @info("Entry", rxni, rxn_id, met_id, stoi)
                
                stoi = nothing
                continue
            end

        end # for w

        stoi = nothing
        sense = -1

    end # for rxni

    unique!(met_ids)

    M = length(met_ids)
    N = length(rxn_ids)
    
    S = spzeros(M, N)
    # S = zeros(M, N)
    for (meti, met_id) in enumerate(met_ids)
        for (rxni, rxn_id) in enumerate(rxn_ids)
            S[meti, rxni] = get(S_dict, (met_id, rxn_id), 0.0)
        end
    end
    
    # bounds
    lb = zeros(N)
    ub = zeros(N) .+ 1000.0

    # checkings
    for _ in 1:3
        for (i, met) in enumerate(met_ids)
            c = count(!iszero, S[i,:])
            if c == 1
                # @warn string("met (", met, ") is involved in only one reaction! It will be deleted")
                S[i,:] .= 0.0
            end
        end
    end
    
    net = MetNet(;S,
        b = zeros(M),
        c =  zeros(N),
        lb, ub,
        rxns = string.(rxn_ids),
        mets = string.(met_ids),
    )
    net = _common_format(net)

    extras!(net, "ATPM", "vATP(n)")
    linear_weights!(net, "vATP(n)", 1.0)

    return net
    
end

function _register_Massucci2013()
    register_network!("Massucci2013", _Massucci2013_builder;
        use_cache = false,
        source = "Massucci, Francesco A., Mauro DiNuzzo, Federico Giove, Bruno Maraviglia, Isaac Perez Castillo, Enzo Marinari, and Andrea De Martino. “Energy Metabolism and Glutamate-Glutamine Cycle in the Brain: A Stoichiometric Modeling Perspective.” BMC Systems Biology 7, no. 1 (October 10, 2013): 103. https://doi.org/10.1186/1752-0509-7-103.", 
        desc = "Network used at Massucci 2013 modeling both neuron and astrocyte metabolism"
    )
end