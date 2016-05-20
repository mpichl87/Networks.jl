__precompile__()
module Networks 
importall Base

const Z0 = 50.0
export Z0

type Network{ T <: Number }
	sparams::Array{ T, 2 }
	measures::Array{ T, 2 }
	labels::Array{ ASCIIString, 1 }
	function Network{ T }( sparams::Array{ T, 2 }, measures::Array{ T, 2 }, labels::Array{ ASCIIString, 1 } )
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
function Network{ T }( sparams::Array{ T, 2 }, measures::Array{ T, 2 }, labels::Array{ ASCIIString, 1 } ) 
	Network{ T }( sparams, measures, labels )
end
function Network{ S, T }( sparams::Array{ S, 2 }, measures::Array{ T, 2 }, labels::Array{ ASCIIString, 1 } )
	sparams, measures = promote( sparams, measures )
	Network{ typeof( sparams ) }( sparams, measures, labels )
end
function Network{ T }( sparams::Array{ T, 2 } )
	Network{ T }( sparams, Array{ T, 2 }(), ASCIIString[] ) 
end
export Network

function ==( nw1::Network, nw2::Network )
	ok = true
	ok &= nw1.sparams == nw2.sparams 
	ok &= ( length( nw1.measures ) == 0 && length( nw2.measures ) == 0 ) || nw1.measures == nw2.measures
	ok &= nw1.labels == nw2.labels
end

function measure( S::Network, label::ASCIIString )
	for ( i, l ) in enumerate( S.labels )
		if l == label
			return S.measures[ i, : ]
		end
	end
end
export measure

voltage( S::Network, component::ASCIIString ) = measure( S, "$( component )-U")
voltage( S::Network ) = component -> voltage( S, component )
export  voltage

current( S::Network, component::ASCIIString ) = measure( S, "$( component )-I")
current( S::Network ) = component -> current( S, component )
export current

power( S::Network, component::ASCIIString ) = voltage( S, component ) * conj( current( S, component ) ) / 2
power( S::Network ) = component -> power( S, component )
export power

function UI( component::ASCIIString = "" )
	UIsp = [  	0.0	1.0;
				1.0	0.0 ]
	UImeas = [ 	√(2Z0) 	√(2Z0);
				√2/√Z0 -√2/√Z0	]
	UIlab = [ "$( component )-U", "$( component )-I" ]
	Network( UIsp, UImeas, UIlab )
end
export UI

function UI2( component::ASCIIString = "" )
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
 
const TEEsp = [	-1.0/3  2.0/3  2.0/3;
				 2.0/3 -1.0/3  2.0/3;
				 2.0/3  2.0/3 -1.0/3 ]
const Tee = Network( TEEsp )
export Tee

shunt( z, z0 ) = Network( reshape( [ ( z - z0 ) / ( z + z0 ) ], 1, 1 ) )
shunt( z ) = shunt( z, Z0 )
export shunt

through( z, z0 ) = Network( 1 / ( z + 2z0 ) * [ z 2z0; 2z0 z ] )
through( z ) = through( z, Z0 )
export through

function connect{ S, T }( SN::Network{ S }, k::Int, TN::Network{ T }, l::Int ) 
	S, T = promote( SN.sparams, TN.sparams )
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
	A, B = promote( SN.measures, TN.measures )
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
	labels 		= [ SN.labels; 	TN.labels ]
	Network( sparams, measures, labels )
end

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

function labels!( nw::Network, prefix::ASCIIString )
	labels = map( l -> "$( prefix )$( l )", nw.labels )
	Network( nw.sparams, nw.measures, labels )
end

instr_shunt(   sh::Network, component::ASCIIString = "" ) = connect( sh, 1, UI( component ), 2 )
instr_through( th::Network, component::ASCIIString = "" ) = connect( connect( th, 1, UI2( component ), 2 ), 1, 3 )
export instr_shunt, instr_through

parallel( sh1::Network, sh2::Network ) = connect( connect( Tee, 1, sh1, 1 ), 1, sh2, 1 )
serial(    th::Network,  sh::Network ) = connect( th, 2, sh, 1 )
export parallel, serial

y( z ) = 1 / z
z( y ) = 1 / y
export y, z

ser_z( z1, z2 ) = z1 + z2
par_y( y1, y2 ) = y1 + y2
ser = ser_z
export ser_z, par_y, ser

par_z( z1, z2 ) = z( par_y( y( z1 ), y( z2 ) ) )
ser_y( y1, y2 ) = y( ser_z( z( y1 ), z( y2 ) ) )
par = par_z
export par_z, ser_y, par

ω( f ) = 2π * f
f( ω ) = ω / 2π
export ω, f

z_l( l, ω ) = ω * l * 1im
y_l( l, ω ) = y( z_l( l, ω ) )
export z_l, y_l

y_c( c, ω ) = ω * c * 1im
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
