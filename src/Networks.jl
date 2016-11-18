__precompile__()
module Networks
importall Base

const Z0 = 50.0
export Z0

type Network{ NumType <: Number }
	sparams::Array{ NumType, 2 }
	measures::Array{ NumType, 2 }
	labels::Array{ String, 1 }
	function Network{ NumType }( sparams::Array{ NumType, 2 }, measures::Array{ NumType, 2 }, labels::Array{ String, 1 } )
		srows = size( sparams, 1 )
		scols = srows == 0 ? 0 : size( sparams, 2 )
		mrows = size( measures, 1 )
		mcols = mrows == 0 ? scols : size( measures, 2 )
		lrows = size( labels, 1 )
		srows != scols ? error( "sparams Matrix not square" ) :
		scols != mcols ? error( "number of columns of sparams and measures differ") :
		mrows != lrows ? error( "number of rows of measures and labels differ") :
		new( sparams, measures, labels )
	end
end
function Network{ NumType }( sparams::Array{ NumType, 2 }, measures::Array{ NumType, 2 }, labels::Array{ String, 1 } )
	Network{ NumType }( sparams, measures, labels )
end
function Network{ NumType1, NumType2 }( sparams::Array{ NumType1, 2 }, measures::Array{ NumType2, 2 }, labels::Array{ String, 1 } )
	sparams, measures = promote( sparams, measures )
	Network{ typeof( sparams ) }( sparams, measures, labels )
end
function Network{ NumType }( sparams::Array{ NumType, 2 } )
	Network{ NumType }( sparams, Array{ NumType, 2 }(), String[] )
end
export Network

"""
		==( nw1, nw2 )

Test equality of two networks.
"""
function ==( nw1::Network, nw2::Network )
	ok = true
	ok &= nw1.sparams == nw2.sparams
	ok &= ( length( nw1.measures ) == 0 && length( nw2.measures ) == 0 ) || nw1.measures == nw2.measures
	ok &= nw1.labels == nw2.labels
end

"""
		measure( S, label )

Extracts labeled measurement from network S.
"""
function measure( S::Network, label::String )
	for ( i, l ) in enumerate( S.labels )
		if l == label
			return S.measures[ i, : ]
		end
	end
end
export measure

"""
		voltage( S, component )

Extracts  voltage measurement labeled with component from network S.
"""
voltage( S::Network, component::String ) = measure( S, "$( component )-U")
voltage( S::Network ) = component -> voltage( S, component )
export  voltage

"""
		current( S, component )

Extracts  current measurement labeled with component from network S.
"""
current( S::Network, component::String ) = measure( S, "$( component )-I")
current( S::Network ) = component -> current( S, component )
export current

"""
		power( S, component )

Computes  power measurement from voltage and current labeled with component from network S.
"""
power( S::Network, component::String ) = voltage( S, component ) * conj( current( S, component ) ) / 2
power( S::Network ) = component -> power( S, component )
export power

"""
		UI( component )

Generates voltage and current measurement ports labeled with component for one port.
"""
function UI( component::String = "" )
	UIsp = [  	0.0	1.0;
				1.0	0.0 ]
	UImeas = [ 	√(2Z0) 	√(2Z0);
				√2/√Z0 -√2/√Z0	]
	UIlab = [ "$( component )-U", "$( component )-I" ]
	Network( UIsp, UImeas, UIlab )
end
export UI

"""
		UI2( component )

Generates voltage and current measurement ports labeled with component for two ports.
"""
function UI2( component::String = "" )
	UI2sp = [ 	0.0			1.0 		0.0 		0.0;
				1.0			0.0			0.0			0.0;
				0.0			0.0			0.0			1.0;
				0.0			0.0			1.0			0.0 ]
	UI2meas = [ 	√(2Z0) 		√(2Z0)		-√(2Z0) 	-√(2Z0);
					1/√(2Z0)	-1/√(2Z0)	1/√(2Z0)	-1/√(2Z0)	]
	UI2lab = [ "$( component )-U", "$( component )-I" ]
	Network( UI2sp, UI2meas, UI2lab )
end
export UI2
"""
		TEEsp

S-parameters for 3-port tee.
"""
const TEEsp = [	-1.0/3  2.0/3  2.0/3;
								 2.0/3 -1.0/3  2.0/3;
				 				 2.0/3  2.0/3 -1.0/3 ]

"""
		TEEsp

3-port tee network.
"""
const Tee = Network( TEEsp )
export Tee

"""
		shunt( z[, z0] )

Generates 1-port shunt network of a shunted impedanze z.
z0 is the system impedance with a default value of 50 Ohms.
"""
shunt( z, z0 ) = Network( reshape( [ ( z - z0 ) / ( z + z0 ) ], 1, 1 ) )
shunt( z ) = shunt( z, Z0 )
export shunt

"""
		through( z[, z0] )

Generates 2-port through network from an impedanze z.
z0 is the system impedance with a default value of 50 Ohms.
"""
through( z, z0 ) = Network( 1 / ( z + 2z0 ) * [ z 2z0; 2z0 z ] )
through( z ) = through( z, Z0 )
export through

"""
		connect( N1, k, N2, l )

Connects port k of network N1 with port l of Network N2 and returns the resulting network.
"""
function connect{ NumType1, NumType2 }( N1::Network{ NumType1 }, k::Int, N2::Network{ NumType2 }, l::Int )
	S, T = promote( N1.sparams, N2.sparams )
	nS = size( S )[ 1 ]
	nT = size( T )[ 1 ]
	nR = nS + nT - 2
	sparams = zeros( eltype( S ), nR, nR )
	Skj = S[ k, : ]
	Tlj = T[ l, : ]
	Sik = S[ :, k]
	Til = T[ :, l]
 	Skk = Skj[ k ]
	Tll = Tlj[ l ]
	invden = 1 / ( 1 - Skk * Tll )
	for i in 1:k - 1
		fak2 = Sik[ i ] * invden
		fak1 =  Tll * fak2
		for j in 1:k - 1
			sparams[ i, 			j			] = S[ i, j ] + Skj[ j ] * fak1
		end
		for j in k + 1:nS
			sparams[ i, 			j - 1 		] = S[ i, j ] + Skj[ j ] * fak1
		end
		for j in 1:l - 1
			sparams[ i, 			j + nS - 1	] = Tlj[ j ] * fak2
		end
		for j in l + 1:nT
			sparams[ i, 			j + nS - 2 	] = Tlj[ j ] * fak2
		end
	end
	for i in k + 1:nS
		fak2 = Sik[ i ] * invden
		fak1 =  Tll * fak2
		for j in 1:k - 1
			sparams[ i - 1, 		j 			] = S[ i, j ] + Skj[ j ] * fak1
		end
			for j in k + 1:nS
			sparams[ i - 1,		j - 1 		] = S[ i, j ] + Skj[ j ] * fak1
		end
		for j in 1:l - 1
			sparams[ i - 1, 		j + nS - 1 	] = Tlj[ j ] * fak2
		end
		for j in l + 1:nT
			sparams[ i - 1, 		j + nS - 2 	] = Tlj[ j ] * fak2
		end
	end
	for i in 1:l - 1
		fak1 = Til[ i ] * invden
		fak2 = Skk * fak1
		for j in 1:k - 1
			sparams[ i + nS - 1, 	j 			] = Skj[ j ]  * fak1
		end
		for j in k + 1:nS
			sparams[ i + nS - 1, 	j - 1 		] = Skj[ j ]  * fak1
		end
		for j in 1:l - 1
			sparams[ i + nS - 1, 	j + nS - 1	] = T[ i, j ] + Tlj[ j ] * fak2
		end
		for j in l + 1:nT
			sparams[ i + nS - 1, 	j + nS - 2 	] = T[ i, j ] + Tlj[ j ] * fak2
		end
	end
	for i in l + 1:nT
		fak1 = Til[ i ] * invden
		fak2 = Skk * fak1
		for j in 1:k - 1
			sparams[ i + nS - 2, 	j 			] = Skj[ j ] * fak1
		end
		for j in k + 1:nS
			sparams[ i + nS - 2, 	j - 1 		] = Skj[ j ] * fak1
		end
		for j in 1:l - 1
			sparams[ i + nS - 2, 	j + nS - 1	] = T[ i, j ] + Tlj[ j ] * fak2
		end
		for j in l + 1:nT
			sparams[ i + nS - 2, 	j + nS - 2 	] = T[ i, j ] + Tlj[ j ] * fak2
		end
	end
	A, B = promote( N1.measures, N2.measures )
	nA = size( A )[ 1 ]
	nB = size( B )[ 1 ]
	nM = nA + nB
	measures = zeros( eltype( A ), nM, nR )
	if nA > 0
		Aik = A[ :, k]
		for i in 1:nA
			fak2 = Aik[ i ] * invden
			fak1 =  Tll * fak2
			for j in 1:k - 1
				measures[ i, 	j 			] = A[ i, j ] + Skj[ j ] * fak1
			end
				for j in k + 1:nS
				measures[ i, 	j - 1 		] = A[ i, j ] + Skj[ j ] * fak1
			end
			for j in 1:l - 1
				measures[ i, 	j + nS - 1 ] = Tlj[ j ] * fak2
			end
			for j in l + 1:nT
				measures[ i, 	j + nS - 2 ] = Tlj[ j ] * fak2
			end
		end
	end
	if nB > 0
		Bil = B[ :, l]
		for i in 1:nB
			fak1 = Bil[ i ] * invden
			fak2 = Skk * fak1
			for j in 1:k - 1
				measures[ i + nA, j 			] = Skj[ j ] * fak1
			end
			for j in k + 1:nS
				measures[ i + nA, j - 1 		] = Skj[ j ] * fak1
			end
			for j in 1:l - 1
				measures[ i + nA, j + nS - 1	] = B[ i, j ] + Tlj[ j ] * fak2
			end
			for j in l + 1:nT
				measures[ i + nA, j + nS - 2 	] = B[ i, j ] + Tlj[ j ] * fak2
			end

		end
	end
	labels 		= [ N1.labels; 	N2.labels ]
	Network( sparams, measures, labels )
end

"""n		connect( SN, k, l )

Connects port k of network SN with port l of the same network and returns the resulting network.
"""
function connect( SN::Network, k::Int, l::Int )
	if k == l
		SN
	else
		if k > l
			l, k = k, l
		end
		S = SN.sparams
		nS = size( S )[ 1 ]
		nR = nS - 2
		sparams = zeros( eltype( S ), nR, nR )
		Skj = S[ k, : ]
		Sik = S[ :, k]
	 	Slj = S[ l, : ]
		Sil = S[ :, l]
		Skl = Skj[ l ]
	 	Skk = Skj[ k ]
		Sll = Slj[ l ]
		Slk = Slj[ k ]
		invden = 1 / ( 1 - Skl ) * ( 1 - Slk ) - Skk * Sll
		for i in 1:k - 1
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				sparams[ i, j 	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				sparams[ i, j - 1 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				sparams[ i, j - 2 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		for i in k + 1:l - 1
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				sparams[ i - 1, j 	] =	( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				sparams[ i - 1, j - 1 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				sparams[ i - 1, j - 2 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		for i in l + 1:nS
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				sparams[ i - 2, j 	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				sparams[ i - 2, j - 1	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				sparams[ i - 2, j - 2	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		A = SN.measures
		nA = size( A )[ 1 ]
		measures = zeros( eltype( A ), nA, nR )
		Aik = A[ :, k]
	 	Ail = A[ :, l]

		for i in 1:nA
			fak1 = Sll * Aik[ i ] + Ail[ i ] * ( 1 - Slk )
			fak2 = Skk * Ail[ i ] + Aik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				measures[ i, j 	] = ( A[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				measures[ i, j - 1	] = ( A[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				measures[ i, j - 2	] = ( A[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		labels 		= SN.labels
		Network( sparams, measures, labels )
	end
end
export connect

"""
		prefix_labels( nw, prefix )

Prefixes all labels of a network with a common prefix.
"""
function prefix_labels!( nw::Network, prefix::String )
	labels = map( l -> "$( prefix )$( l )", nw.labels )
	Network( nw.sparams, nw.measures, labels )
end

"""
		instr_shunt( sh, component )

Instrumentates a shunt 1-port network. Adds voltage and current measurement ports labeled with component.
"""
instr_shunt(   sh::Network, component::String = "" ) = connect( sh, 1, UI( component ), 2 )

"""
		instr_through( sh, component )

Instrumentates a through 2-port network. Adds voltage and current measurement ports labeled with component.
"""
instr_through( th::Network, component::String = "" ) = connect( connect( th, 1, UI2( component ), 2 ), 1, 3 )
export instr_shunt, instr_through

"""
		parallel( sh1, sh2 )

Connects two shunt networks in parallel.
"""
parallel( sh1::Network, sh2::Network ) = connect( connect( Tee, 1, sh1, 1 ), 1, sh2, 1 )

"""
		serial( th, sh )

Connects a through network in series to a shunt network.
"""
serial(    th::Network,  sh::Network ) = connect( th, 2, sh, 1 )
export parallel, serial

"""
		y( z )

Calculates the admitance from an impedance.
"""
y( z ) = 1 / z

"""
		z( y )

Calculates the impedance from an admitance.
"""
z( y ) = 1 / y
export y, z

"""
		ser_z( z1, z2 )

Calculates the impedance of a serial connection of two impedances.
"""
ser_z( z1, z2 ) = z1 + z2

"""
		par_y( y1, y2 )

Calculates the admitance of a parallel connection of two admitances.
"""
par_y( y1, y2 ) = y1 + y2
ser = ser_z
export ser_z, par_y, ser

"""
		par_z( z1, z2 )

Calculates the impedance of a parallel connection of two impedances.
"""
par_z( z1, z2 ) = z( par_y( y( z1 ), y( z2 ) ) )

"""
		ser_y( y1, y2 )

Calculates the admitance of a serial connection of two admitances.
"""
ser_y( y1, y2 ) = y( ser_z( z( y1 ), z( y2 ) ) )
par = par_z
export par_z, ser_y, par

"""
		ω( f )

Calculates the angular frequency from frequency.
"""
ω( f ) = 2π * f

"""
		f( ω )

Calculates the frequency from angular frequency.
"""
f( ω ) = ω / 2π
export ω, f

"""
		z_l( l, ω )

Calculates the impedance of an ideal inductor with inductance l at the  angular frequency ω.
"""
z_l( l, ω ) = ω * l * 1im


"""
		y_l( l, ω )

Calculates the admitance of an ideal inductor with inductance l at the  angular frequency ω.
"""
y_l( l, ω ) = y( z_l( l, ω ) )
export z_l, y_l

"""
		y_c( c, ω )

Calculates the admitance of an ideal capacitor with capacitance c at the  angular frequency ω.
"""
y_c( c, ω ) = ω * c * 1im

"""
		z_c( c, ω )

Calculates the impedance of an ideal capacitor with capacitance c at the  angular frequency ω.
"""
z_c( c, ω ) = z( y_c( c, ω ) )
export y_c, z_c

rcl( r, c, l, ω ) =
	ser(
		par( r, z_c( c, ω ) ),
		z_l( l, ω )
	)

rcl( r, c, l ) = ω -> rcl( r, c, l, ω )
export rcl


cap( c, r, ω ) = ser( r, z_c( c, ω ) )
ind( l, r, ω ) = ser( r, z_l( l, ω ) )
export cap, ind

function pin( Rsf, Ct )
	return function ( state, ω )
		state ? Rsf : z_c( Ct, ω )
	end
end
export pin

end
