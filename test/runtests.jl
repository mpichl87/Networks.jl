using Base.Test
using Networks 

@test shunt( Z0 ) == Network( reshape( [ 0.0 ], 1, 1 ) )
# @test shunt( 0 )  = Network( [ -1 ] )

# @test through( 0 ) = Network( [  0 1;
# 						1 0 ] )

# @test Networks.connect( through( 0 ), 2, shunt( Z0 ), 1 ) = Network( shunt( Z0 ) )

