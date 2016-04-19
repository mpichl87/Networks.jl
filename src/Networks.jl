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

# http://qucs.sourceforge.net/docs/technical/technical.pdf p9 (1.14) (1.15)
function connect( S, k, T, l ) 
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
# http://qucs.sourceforge.net/docs/technical/technical.pdf p10 (1.16)
function connect( S, k, l ) 
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
export connect

instr_shunt( sh ) 	= connect( sh, 1, UI, 2 )
instr_through( th ) = connect( connect( th, 1, UI2, 2 ), 1, 3 )
export instr_shunt, instr_through

parallel( sh1, sh2 ) = nw.connect( nw.connect( Tee, 1, sh1, 1 ), 1, sh2, 1 )
serial( th, sh )   	 = nw.connect( th, 2, sh, 1 )
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
