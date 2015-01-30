pro make_idl_aacgm64

  cd,'.',current=dir

  files=['idl_aacgm','idl_aacgm_call','aacgm','altitude_to_cgm','cgm_to_altitude','coeff',$
    'convert_geo_coord','convert_geo_vec','eqn_of_time','math','mlt','mlt1','rylm','solar_loc']

  ld='gcc -shared -o %L %O %X'

  make_dll,files,'IDL_Load',compile_directory=dir+'/src',input_directory=dir+'/src', $
    output_directory=dir,ld=ld

return
end
