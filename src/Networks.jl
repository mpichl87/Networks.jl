__precompile__()
module Networks 

typealias F Float64
typealias C Complex{ F }

const Z0 = 50.0
export Z0

const UI = [ 	0.0		1.0 	0.0 0.0;
				1.0		0.0		0.0	0.0;
				√(2Z0) 	√(2Z0)	0.0	0.0;
				√2/√Z0 -√2/√Z0	0.0	0.0	]
export UI

const UI2 = [ 	0.0			1.0 		0.0 		0.0 		0.0	0.0;
				1.0			0.0			0.0			0.0 		0.0	0.0;
				0.0			0.0			0.0			1.0 		0.0	0.0;
				0.0			0.0			1.0			0.0 		0.0	0.0;
				√(2Z0) 		√(2Z0)		-√(2Z0) 	-√(2Z0) 	0.0	0.0;
				1/√(2Z0)	-1/√(2Z0)	1/√(2Z0)	-1/√(2Z0) 	0.0	0.0  ]
export UI2

 
const Tee = [ 	-1.0/3  2.0/3  2.0/3;
				 2.0/3 -1.0/3  2.0/3;
				 2.0/3  2.0/3 -1.0/3 ]
export Tee

shunt( z, z0 ) = [ ( z - z0 ) / ( z + z0 ) ]
shunt( z ) = shunt( z, Z0 )
export shunt

through( z, z0 ) = 1 / ( z + 2z0 ) * [ z 2z0; 2z0 z ]
through( z ) = through( z, Z0 )
export through

type Network
	sparams
	measures
	labels
end

Network( sparams::Array ) = Network( sparams, [], [] )
export Network

const UI2nw = [0.0			1.0 		0.0 		0.0;
				1.0			0.0			0.0			0.0;
				0.0			0.0			0.0			1.0;
				0.0			0.0			1.0			0.0 ]
const UI2meas = [ 	√(2Z0) 		√(2Z0)		-√(2Z0) 	-√(2Z0);  
					1/√(2Z0)	-1/√(2Z0)	1/√(2Z0)	-1/√(2Z0)	]
const UI2lab = [ "U", "I" ]

const UI2NW = Network( UI2nw, UI2meas, UI2lab )
export UI2NW



# http://qucs.sourceforge.net/docs/technical/technical.pdf p9 (1.14) (1.15)
function connect( S::Array, k::Int, T::Array, l::Int ) 
	nS = size( S )[ 1 ]
	nT = size( T )[ 1 ]
	nR = nS + nT - 2
	R = zeros( C, nR, nR )
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
			R[ i, 			j			] = S[ i, j ] + Skj[ j ] * fak1
		end
		for j in k + 1:nS
			R[ i, 			j - 1 		] = S[ i, j ] + Skj[ j ] * fak1
		end
		for j in 1:l - 1
			R[ i, 			j + nS - 1	] = Tlj[ j ] * fak2
		end
		for j in l + 1:nT 
			R[ i, 			j + nS - 2 	] = Tlj[ j ] * fak2
		end
	end
	for i in k + 1:nS
		fak2 = Sik[ i ] * invden
		fak1 =  Tll * fak2
		for j in 1:k - 1
			R[ i - 1, 		j 			] = S[ i, j ] + Skj[ j ] * fak1
		end
			for j in k + 1:nS
			R[ i - 1,		j - 1 		] = S[ i, j ] + Skj[ j ] * fak1
		end
		for j in 1:l - 1
			R[ i - 1, 		j + nS - 1 	] = Tlj[ j ] * fak2
		end
		for j in l + 1:nT 
			R[ i - 1, 		j + nS - 2 	] = Tlj[ j ] * fak2
		end
	end
	for i in 1:l - 1
		fak1 = Til[ i ] * invden 
		fak2 = Skk * fak1
		for j in 1:k - 1
			R[ i + nS - 1, 	j 			] = Skj[ j ]  * fak1
		end
		for j in k + 1:nS
			R[ i + nS - 1, 	j - 1 		] = Skj[ j ]  * fak1
		end
		for j in 1:l - 1
			R[ i + nS - 1, 	j + nS - 1	] = T[ i, j ] + Tlj[ j ] * fak2
		end
		for j in l + 1:nT
			R[ i + nS - 1, 	j + nS - 2 	] = T[ i, j ] + Tlj[ j ] * fak2
		end
		
	end
	for i in l + 1:nT
		fak1 = Til[ i ] * invden 
		fak2 = Skk * fak1
		for j in 1:k - 1
			R[ i + nS - 2, 	j 			] = Skj[ j ] * fak1
		end
		for j in k + 1:nS
			R[ i + nS - 2, 	j - 1 		] = Skj[ j ] * fak1
		end
		for j in 1:l - 1
			R[ i + nS - 2, 	j + nS - 1	] = T[ i, j ] + Tlj[ j ] * fak2
		end
		for j in l + 1:nT
			R[ i + nS - 2, 	j + nS - 2 	] = T[ i, j ] + Tlj[ j ] * fak2
		end
		
	end
	R
end

function connect( SN::Network, k::Int, TN::Network, l::Int ) 
	S = SN.sparams
	T = TN.sparams
	nS = size( S )[ 1 ]
	nT = size( T )[ 1 ]
	nR = nS + nT - 2
	sparams = zeros( C, nR, nR )
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
	A = SN.measures
	B = TN.measures
	nA = size( A )[ 1 ]
	nB = size( B )[ 1 ]
	nM = nA + nB
	measures = zeros( C, nM, nR )
	Aik = A[ :, k]
	Bil = B[ :, l]
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

	measures 	= R
	labels 		= [ S.labels; 	T.labels ]
	Network( sparams, measures, labels )
end


# http://qucs.sourceforge.net/docs/technical/technical.pdf p10 (1.16)
function connect( S::Array, k::Int, l::Int ) 
	if k == l
		S
	else
		if k > l
			l, k = k, l
		end
		nS = size( S )[ 1 ]
		nR = nS - 2
		R = zeros( C, nR, nR )
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
				R[ i, j 	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i, j - 1 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i, j - 2 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		for i in k + 1:l - 1
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				R[ i - 1, j 	] =	( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i - 1, j - 1 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i - 1, j - 2 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		for i in l + 1:nS
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				R[ i - 2, j 	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i - 2, j - 1	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i - 2, j - 2	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		R
	end
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
		R = zeros( C, nR, nR )
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
				R[ i, j 	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i, j - 1 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i, j - 2 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		for i in k + 1:l - 1
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				R[ i - 1, j 	] =	( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i - 1, j - 1 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i - 1, j - 2 ] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		for i in l + 1:nS
			fak1 = Sll * Sik[ i ] + Sil[ i ] * ( 1 - Slk )
			fak2 = Skk * Sil[ i ] + Sik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				R[ i - 2, j 	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i - 2, j - 1	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i - 2, j - 2	] = ( S[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		sparams = R
		A = SN.measures
		nA = size( A )[ 1 ]
		R = zeros( C, nA, nR )
		Aik = A[ :, k]
	 	Ail = A[ :, l]
		
		for i in 1:nA
			fak1 = Sll * Aik[ i ] + Ail[ i ] * ( 1 - Slk )
			fak2 = Skk * Ail[ i ] + Aik[ i ] * ( 1 - Skl )
			for j in 1:k - 1
				R[ i, j 	] = ( A[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in k + 1:l - 1
				R[ i, j - 1	] = ( A[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
			for j in l + 1:nS
				R[ i, j - 2	] = ( A[ i, j ] + Skj[ j ] * fak1 +	Slj[ j ] * fak2 ) * invden
			end
		end
		
		measures 	= R
		labels 		= SN.labels
		Network( sparams, measures, labels )
	end
end

export connect

instr_shunt( sh ) 	= connect( sh, 1, UI, 2 )
instr_through( th ) = connect( connect( th, 1, UI2, 2 ), 1, 3 )
export instr_shunt, instr_through

parallel( sh1, sh2 ) = connect( connect( Tee, 1, sh1, 1 ), 1, sh2, 1 )
serial( th, sh )   	 = connect( th, 2, sh, 1 )
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

end
