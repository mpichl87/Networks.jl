using Base.Test
using Networks

@test shunt( Z0 ) 	== Network( reshape( [ 0.0 ], 1, 1 ) )
@test shunt( 0.0 )  == Network( reshape( [ -1.0 ], 1, 1 ) )

@test through( 0 ) 	== Network( [ 	0 1;
									1 0 ] )

@test Networks.conn( through( 0.0 ), 2, shunt( Z0 ), 1 ) == Network( shunt( Z0 ) )


@test_throws ErrorException Network( [ 1 0; 0 1; 0 1 ] )

@test Network( [ 1 0; 0 1 ] ).sparams == [ 1 0; 0 1 ]

info( "pwd: $( pwd() )" )
info( "project dir: $( readdir( Pkg.dir( "Networks" ) ) )" )
info( "docs dir: $( readdir( "$(Pkg.dir( "Networks" ))/docs" ) )" )
