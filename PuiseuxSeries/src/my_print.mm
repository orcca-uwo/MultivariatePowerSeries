export my_print := proc( a :: seq( name ), $ )
  local L := [ a ];
  local u;
  print( seq( u = eval(u,2), u = L ) );
end proc;

export my_print2 := proc( a :: uneval, b :: uneval, c :: uneval, d :: uneval, e :: uneval, $ )
  local L, u;
  L := [ _passed ];
  print( seq( u = eval(u,2), u = L ) );
end proc;