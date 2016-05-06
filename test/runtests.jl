using Base.Test
using Networks 

@test shunt( Z0 ) ≈ [ 0 ]
@test shunt( 0 )  ≈ [ -1 ]

@test through( 0 ) ≈ [  0 1;
						1 0 ]

@test Networks.connect( through( 0 ), 2, shunt( Z0 ), 1 ) ≈ shunt( Z0 )
