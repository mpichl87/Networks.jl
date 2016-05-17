using Base.Test
using Networks 

@test shunt( Z0 ) 	== Network( reshape( [ 0.0 ], 1, 1 ) )
@test shunt( 0.0 )  == Network( reshape( [ -1.0 ], 1, 1 ) )

@test through( 0 ) 	== Network( [ 	0 1;
									1 0 ] )

@test Networks.connect( through( 0.0 ), 2, shunt( Z0 ), 1 ) == Network( shunt( Z0 ) )


@test_throws ErrorException Network( [ 1 0; 0 1; 0 1 ] )

@test Network( [ 1 0; 0 1 ] ).sparams == [ 1 0; 0 1 ]

# Network( [ 0.0 1.0;
# 				 1.0 0.0 ],
# 				 reshape( [ 0.0 ], 1, 1 ), 
# 				 [ "test" ] ) == 1

