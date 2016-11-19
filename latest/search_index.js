var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Networks.jl Documentation",
    "title": "Networks.jl Documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Networks.jl-Documentation-1",
    "page": "Networks.jl Documentation",
    "title": "Networks.jl Documentation",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#Status-1",
    "page": "Networks.jl Documentation",
    "title": "Status",
    "category": "section",
    "text": "(Image: Travis Status) (Image: Coverage Status)"
},

{
    "location": "index.html#Description-1",
    "page": "Networks.jl Documentation",
    "title": "Description",
    "category": "section",
    "text": "One method for S-parameter simulation of linear networks is to connect the ports of smaller S-matrices. This method is used in qucs.There are two cases to consider:Connect two ports of the same network\nIf two ports of the same circuit S are connected, the new S-parameters are\nS^prime_ij = S_ij + frac S_kj  S_il  ( 1  S_lk ) + S_lj  S_ik  ( 1  S_kl ) + S_kj  S_ll  S_ik + S_lj  S_kk  S_il   ( 1  S_kl )  ( 1  S_lk )  S_kk  S_ll \nConnect two ports of two different networks:\nConnecting port k of circuit S with port l of circuit T, the new S-parameters are\nS^prime_ij = S_ij +  frac S_kj  T_ll  S_ik  1  S_kk  T_ll \nwith i and j both being ports of S. Furthermore, it is\nS^prime_mj = frac S_kj  T_ml  1  S_kk  T_ll \nwith m being a port of the circuit T.S-parameters are \"Black Boxes\". No internal state is visible. To make measurement inside a circuit, the matrices are augmented with ports to output internal signals ( voltage, current ). Because there is no input from this ports, the resulting matrix gets long but narrow. This output-only part is held in a different matrix as the normal s-parameters. Additionally, the lines of this matrix can be labeled with a name, e.g the name of the node + \"_V\" for voltage.To simulate a network and measure some values:Construct the networks for the individual components.\nAugment the matrix for the ports you want to measure with output-only ports.\nConnect the networks together, until only one port for each source remains.This results in a network with a small s-parameter matrix, which describes the flow of power between all sources through the circuit. The additional lines in the measurement matrix describe the measured signals in dependence of the power input from each (real) port."
},

{
    "location": "index.html#Networks.Tee",
    "page": "Networks.jl Documentation",
    "title": "Networks.Tee",
    "category": "Constant",
    "text": "	Tee\n\n3-port tee network.\n\n\n\n"
},

{
    "location": "index.html#Networks.UI",
    "page": "Networks.jl Documentation",
    "title": "Networks.UI",
    "category": "Function",
    "text": "	UI( component )\n\nGenerates voltage and current measurement ports labeled with component for one port.\n\n\n\n"
},

{
    "location": "index.html#Networks.UI2",
    "page": "Networks.jl Documentation",
    "title": "Networks.UI2",
    "category": "Function",
    "text": "	UI2( component )\n\nGenerates voltage and current measurement ports labeled with component for two ports.\n\n\n\n"
},

{
    "location": "index.html#Networks.conn-Tuple{Networks.Network,Int64,Int64}",
    "page": "Networks.jl Documentation",
    "title": "Networks.conn",
    "category": "Method",
    "text": "	conn( SN, k, l )\n\nConnects port k of network SN with port l of the same network and returns the resulting network.\n\n\n\n"
},

{
    "location": "index.html#Networks.conn-Tuple{Networks.Network{NumType1},Int64,Networks.Network{NumType2},Int64}",
    "page": "Networks.jl Documentation",
    "title": "Networks.conn",
    "category": "Method",
    "text": "	conn( N1, k, N2, l )\n\nConnects port k of network N1 with port l of Network N2 and returns the resulting network.\n\n\n\n"
},

{
    "location": "index.html#Networks.current-Tuple{Networks.Network,String}",
    "page": "Networks.jl Documentation",
    "title": "Networks.current",
    "category": "Method",
    "text": "	current( S, component )\n\nExtracts  current measurement labeled with component from network S.\n\n\n\n"
},

{
    "location": "index.html#Networks.f-Tuple{Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.f",
    "category": "Method",
    "text": "	f( ω )\n\nCalculates the frequency from angular frequency.\n\n\n\n"
},

{
    "location": "index.html#Networks.instr_shunt",
    "page": "Networks.jl Documentation",
    "title": "Networks.instr_shunt",
    "category": "Function",
    "text": "	instr_shunt( sh, component )\n\nInstrumentates a shunt 1-port network. Adds voltage and current measurement ports labeled with component.\n\n\n\n"
},

{
    "location": "index.html#Networks.instr_through",
    "page": "Networks.jl Documentation",
    "title": "Networks.instr_through",
    "category": "Function",
    "text": "	instr_through( sh, component )\n\nInstrumentates a through 2-port network. Adds voltage and current measurement ports labeled with component.\n\n\n\n"
},

{
    "location": "index.html#Networks.measure-Tuple{Networks.Network,String}",
    "page": "Networks.jl Documentation",
    "title": "Networks.measure",
    "category": "Method",
    "text": "	measure( S, label )\n\nExtracts labeled measurement from network S.\n\n\n\n"
},

{
    "location": "index.html#Networks.par_y-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.par_y",
    "category": "Method",
    "text": "	par_y( y1, y2 )\n\nCalculates the admitance of a parallel connection of two admitances.\n\n\n\n"
},

{
    "location": "index.html#Networks.par_z-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.par_z",
    "category": "Method",
    "text": "	par_z( z1, z2 )\n\nCalculates the impedance of a parallel connection of two impedances.\n\n\n\n"
},

{
    "location": "index.html#Networks.parallel-Tuple{Networks.Network,Networks.Network}",
    "page": "Networks.jl Documentation",
    "title": "Networks.parallel",
    "category": "Method",
    "text": "	parallel( sh1, sh2 )\n\nConnects two shunt networks in parallel.\n\n\n\n"
},

{
    "location": "index.html#Networks.power-Tuple{Networks.Network,String}",
    "page": "Networks.jl Documentation",
    "title": "Networks.power",
    "category": "Method",
    "text": "	power( S, component )\n\nComputes  power measurement from voltage and current labeled with component from network S.\n\n\n\n"
},

{
    "location": "index.html#Networks.ser_y-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.ser_y",
    "category": "Method",
    "text": "	ser_y( y1, y2 )\n\nCalculates the admitance of a serial connection of two admitances.\n\n\n\n"
},

{
    "location": "index.html#Networks.ser_z-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.ser_z",
    "category": "Method",
    "text": "	ser_z( z1, z2 )\n\nCalculates the impedance of a serial connection of two impedances.\n\n\n\n"
},

{
    "location": "index.html#Networks.serial-Tuple{Networks.Network,Networks.Network}",
    "page": "Networks.jl Documentation",
    "title": "Networks.serial",
    "category": "Method",
    "text": "	serial( th, sh )\n\nConnects a through network in series to a shunt network.\n\n\n\n"
},

{
    "location": "index.html#Networks.shunt-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.shunt",
    "category": "Method",
    "text": "	shunt( z[, z0] )\n\nGenerates 1-port shunt network of a shunted impedanze z. z0 is the system impedance with a default value of 50 Ohms.\n\n\n\n"
},

{
    "location": "index.html#Networks.through-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.through",
    "category": "Method",
    "text": "	through( z[, z0] )\n\nGenerates 2-port through network from an impedanze z. z0 is the system impedance with a default value of 50 Ohms.\n\n\n\n"
},

{
    "location": "index.html#Networks.voltage-Tuple{Networks.Network,String}",
    "page": "Networks.jl Documentation",
    "title": "Networks.voltage",
    "category": "Method",
    "text": "	voltage( S, component )\n\nExtracts  voltage measurement labeled with component from network S.\n\n\n\n"
},

{
    "location": "index.html#Networks.y-Tuple{Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.y",
    "category": "Method",
    "text": "	y( z )\n\nCalculates the admitance from an impedance.\n\n\n\n"
},

{
    "location": "index.html#Networks.y_c-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.y_c",
    "category": "Method",
    "text": "	y_c( c, ω )\n\nCalculates the admitance of an ideal capacitor with capacitance c at the  angular frequency ω.\n\n\n\n"
},

{
    "location": "index.html#Networks.y_l-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.y_l",
    "category": "Method",
    "text": "	y_l( l, ω )\n\nCalculates the admitance of an ideal inductor with inductance l at the  angular frequency ω.\n\n\n\n"
},

{
    "location": "index.html#Networks.z-Tuple{Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.z",
    "category": "Method",
    "text": "	z( y )\n\nCalculates the impedance from an admitance.\n\n\n\n"
},

{
    "location": "index.html#Networks.z_c-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.z_c",
    "category": "Method",
    "text": "	z_c( c, ω )\n\nCalculates the impedance of an ideal capacitor with capacitance c at the  angular frequency ω.\n\n\n\n"
},

{
    "location": "index.html#Networks.z_l-Tuple{Any,Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.z_l",
    "category": "Method",
    "text": "	z_l( l, ω )\n\nCalculates the impedance of an ideal inductor with inductance l at the  angular frequency ω.\n\n\n\n"
},

{
    "location": "index.html#Networks.ω-Tuple{Any}",
    "page": "Networks.jl Documentation",
    "title": "Networks.ω",
    "category": "Method",
    "text": "	ω( f )\n\nCalculates the angular frequency from frequency.\n\n\n\n"
},

{
    "location": "index.html#Networks.TEEsp",
    "page": "Networks.jl Documentation",
    "title": "Networks.TEEsp",
    "category": "Constant",
    "text": "	TEEsp\n\nS-parameters for 3-port tee.\n\n\n\n"
},

{
    "location": "index.html#Base.:==-Tuple{Networks.Network,Networks.Network}",
    "page": "Networks.jl Documentation",
    "title": "Base.:==",
    "category": "Method",
    "text": "	==( nw1, nw2 )\n\nTest equality of two networks.\n\n\n\n"
},

{
    "location": "index.html#Networks.prefix_labels!-Tuple{Networks.Network,String}",
    "page": "Networks.jl Documentation",
    "title": "Networks.prefix_labels!",
    "category": "Method",
    "text": "	prefix_labels!( nw, prefix )\n\nPrefixes all labels of a network with a common prefix.\n\n\n\n"
},

{
    "location": "index.html#API-doc-1",
    "page": "Networks.jl Documentation",
    "title": "API doc",
    "category": "section",
    "text": "Modules = [Networks]"
},

{
    "location": "index.html#Index-1",
    "page": "Networks.jl Documentation",
    "title": "Index",
    "category": "section",
    "text": ""
},

]}
