      program fgennc
      include 'netcdf.inc'
* error status return
      integer  iret
* netCDF id
      integer  ncid
* dimension ids
      integer  Dr_dim
      integer  D1_dim
      integer  D2_dim
      integer  D3_dim
      integer  dim_dash_name_dash_dashes_dim
      integer  dim_dot_name_dot_dots_dim
* dimension lengths
      integer  Dr_len
      integer  D1_len
      integer  D2_len
      integer  D3_len
      integer  dim_dash_name_dash_dashes_len
      integer  dim_dot_name_dot_dots_len
      parameter (Dr_len = NF_UNLIMITED)
      parameter (D1_len = 1)
      parameter (D2_len = 2)
      parameter (D3_len = 3)
      parameter (dim_dash_name_dash_dashes_len = 4)
      parameter (dim_dot_name_dot_dots_len = 5)
* variable ids
      integer  c_id
      integer  b_id
      integer  s_id
      integer  i_id
      integer  f_id
      integer  d_id
      integer  cr_id
      integer  br_id
      integer  sr_id
      integer  ir_id
      integer  fr_id
      integer  dr_id
      integer  c1_id
      integer  b1_id
      integer  s1_id
      integer  i1_id
      integer  f1_id
      integer  d1_id
      integer  c2_id
      integer  b2_id
      integer  s2_id
      integer  i2_id
      integer  f2_id
      integer  d2_id
      integer  c3_id
      integer  b3_id
      integer  s3_id
      integer  i3_id
      integer  f3_id
      integer  d3_id
      integer  cr1_id
      integer  br2_id
      integer  sr3_id
      integer  f11_id
      integer  d12_id
      integer  c13_id
      integer  s21_id
      integer  i22_id
      integer  f23_id
      integer  c31_id
      integer  b32_id
      integer  s33_id
      integer  sr11_id
      integer  ir12_id
      integer  fr13_id
      integer  cr21_id
      integer  br22_id
      integer  sr23_id
      integer  fr31_id
      integer  dr32_id
      integer  cr33_id
      integer  c111_id
      integer  b112_id
      integer  s113_id
      integer  f121_id
      integer  d122_id
      integer  c123_id
      integer  s131_id
      integer  i132_id
      integer  f133_id
      integer  f211_id
      integer  d212_id
      integer  c213_id
      integer  s221_id
      integer  i222_id
      integer  f223_id
      integer  c231_id
      integer  b232_id
      integer  s233_id
      integer  s311_id
      integer  i312_id
      integer  f313_id
      integer  var_dash_name_dash_dashes_id
      integer  var_dot_name_dot_dots_id
* rank (number of dimensions) for each variable
      integer  c_rank
      integer  b_rank
      integer  s_rank
      integer  i_rank
      integer  f_rank
      integer  d_rank
      integer  cr_rank
      integer  br_rank
      integer  sr_rank
      integer  ir_rank
      integer  fr_rank
      integer  dr_rank
      integer  c1_rank
      integer  b1_rank
      integer  s1_rank
      integer  i1_rank
      integer  f1_rank
      integer  d1_rank
      integer  c2_rank
      integer  b2_rank
      integer  s2_rank
      integer  i2_rank
      integer  f2_rank
      integer  d2_rank
      integer  c3_rank
      integer  b3_rank
      integer  s3_rank
      integer  i3_rank
      integer  f3_rank
      integer  d3_rank
      integer  cr1_rank
      integer  br2_rank
      integer  sr3_rank
      integer  f11_rank
      integer  d12_rank
      integer  c13_rank
      integer  s21_rank
      integer  i22_rank
      integer  f23_rank
      integer  c31_rank
      integer  b32_rank
      integer  s33_rank
      integer  sr11_rank
      integer  ir12_rank
      integer  fr13_rank
      integer  cr21_rank
      integer  br22_rank
      integer  sr23_rank
      integer  fr31_rank
      integer  dr32_rank
      integer  cr33_rank
      integer  c111_rank
      integer  b112_rank
      integer  s113_rank
      integer  f121_rank
      integer  d122_rank
      integer  c123_rank
      integer  s131_rank
      integer  i132_rank
      integer  f133_rank
      integer  f211_rank
      integer  d212_rank
      integer  c213_rank
      integer  s221_rank
      integer  i222_rank
      integer  f223_rank
      integer  c231_rank
      integer  b232_rank
      integer  s233_rank
      integer  s311_rank
      integer  i312_rank
      integer  f313_rank
      integer  var_dash_name_dash_dashes_rank
      integer  var_dot_name_dot_dots_rank
      parameter (c_rank = 0)
      parameter (b_rank = 0)
      parameter (s_rank = 0)
      parameter (i_rank = 0)
      parameter (f_rank = 0)
      parameter (d_rank = 0)
      parameter (cr_rank = 1)
      parameter (br_rank = 1)
      parameter (sr_rank = 1)
      parameter (ir_rank = 1)
      parameter (fr_rank = 1)
      parameter (dr_rank = 1)
      parameter (c1_rank = 1)
      parameter (b1_rank = 1)
      parameter (s1_rank = 1)
      parameter (i1_rank = 1)
      parameter (f1_rank = 1)
      parameter (d1_rank = 1)
      parameter (c2_rank = 1)
      parameter (b2_rank = 1)
      parameter (s2_rank = 1)
      parameter (i2_rank = 1)
      parameter (f2_rank = 1)
      parameter (d2_rank = 1)
      parameter (c3_rank = 1)
      parameter (b3_rank = 1)
      parameter (s3_rank = 1)
      parameter (i3_rank = 1)
      parameter (f3_rank = 1)
      parameter (d3_rank = 1)
      parameter (cr1_rank = 2)
      parameter (br2_rank = 2)
      parameter (sr3_rank = 2)
      parameter (f11_rank = 2)
      parameter (d12_rank = 2)
      parameter (c13_rank = 2)
      parameter (s21_rank = 2)
      parameter (i22_rank = 2)
      parameter (f23_rank = 2)
      parameter (c31_rank = 2)
      parameter (b32_rank = 2)
      parameter (s33_rank = 2)
      parameter (sr11_rank = 3)
      parameter (ir12_rank = 3)
      parameter (fr13_rank = 3)
      parameter (cr21_rank = 3)
      parameter (br22_rank = 3)
      parameter (sr23_rank = 3)
      parameter (fr31_rank = 3)
      parameter (dr32_rank = 3)
      parameter (cr33_rank = 3)
      parameter (c111_rank = 3)
      parameter (b112_rank = 3)
      parameter (s113_rank = 3)
      parameter (f121_rank = 3)
      parameter (d122_rank = 3)
      parameter (c123_rank = 3)
      parameter (s131_rank = 3)
      parameter (i132_rank = 3)
      parameter (f133_rank = 3)
      parameter (f211_rank = 3)
      parameter (d212_rank = 3)
      parameter (c213_rank = 3)
      parameter (s221_rank = 3)
      parameter (i222_rank = 3)
      parameter (f223_rank = 3)
      parameter (c231_rank = 3)
      parameter (b232_rank = 3)
      parameter (s233_rank = 3)
      parameter (s311_rank = 3)
      parameter (i312_rank = 3)
      parameter (f313_rank = 3)
      parameter (var_dash_name_dash_dashes_rank = 0)
      parameter (var_dot_name_dot_dots_rank = 0)
* variable shapes
      integer  cr_dims(cr_rank)
      integer  br_dims(br_rank)
      integer  sr_dims(sr_rank)
      integer  ir_dims(ir_rank)
      integer  fr_dims(fr_rank)
      integer  dr_dims(dr_rank)
      integer  c1_dims(c1_rank)
      integer  b1_dims(b1_rank)
      integer  s1_dims(s1_rank)
      integer  i1_dims(i1_rank)
      integer  f1_dims(f1_rank)
      integer  d1_dims(d1_rank)
      integer  c2_dims(c2_rank)
      integer  b2_dims(b2_rank)
      integer  s2_dims(s2_rank)
      integer  i2_dims(i2_rank)
      integer  f2_dims(f2_rank)
      integer  d2_dims(d2_rank)
      integer  c3_dims(c3_rank)
      integer  b3_dims(b3_rank)
      integer  s3_dims(s3_rank)
      integer  i3_dims(i3_rank)
      integer  f3_dims(f3_rank)
      integer  d3_dims(d3_rank)
      integer  cr1_dims(cr1_rank)
      integer  br2_dims(br2_rank)
      integer  sr3_dims(sr3_rank)
      integer  f11_dims(f11_rank)
      integer  d12_dims(d12_rank)
      integer  c13_dims(c13_rank)
      integer  s21_dims(s21_rank)
      integer  i22_dims(i22_rank)
      integer  f23_dims(f23_rank)
      integer  c31_dims(c31_rank)
      integer  b32_dims(b32_rank)
      integer  s33_dims(s33_rank)
      integer  sr11_dims(sr11_rank)
      integer  ir12_dims(ir12_rank)
      integer  fr13_dims(fr13_rank)
      integer  cr21_dims(cr21_rank)
      integer  br22_dims(br22_rank)
      integer  sr23_dims(sr23_rank)
      integer  fr31_dims(fr31_rank)
      integer  dr32_dims(dr32_rank)
      integer  cr33_dims(cr33_rank)
      integer  c111_dims(c111_rank)
      integer  b112_dims(b112_rank)
      integer  s113_dims(s113_rank)
      integer  f121_dims(f121_rank)
      integer  d122_dims(d122_rank)
      integer  c123_dims(c123_rank)
      integer  s131_dims(s131_rank)
      integer  i132_dims(i132_rank)
      integer  f133_dims(f133_rank)
      integer  f211_dims(f211_rank)
      integer  d212_dims(d212_rank)
      integer  c213_dims(c213_rank)
      integer  s221_dims(s221_rank)
      integer  i222_dims(i222_rank)
      integer  f223_dims(f223_rank)
      integer  c231_dims(c231_rank)
      integer  b232_dims(b232_rank)
      integer  s233_dims(s233_rank)
      integer  s311_dims(s311_rank)
      integer  i312_dims(i312_rank)
      integer  f313_dims(f313_rank)
* data variables
      integer  b
      integer  s
      integer  i
      real  f
      double precision  d
      integer  b1(D1_len)
      integer  s1(D1_len)
      integer  i1(D1_len)
      real  f1(D1_len)
      double precision  d1(D1_len)
      integer  b2(D2_len)
      integer  s2(D2_len)
      integer  i2(D2_len)
      real  f2(D2_len)
      double precision  d2(D2_len)
      integer  b3(D3_len)
      integer  s3(D3_len)
      integer  i3(D3_len)
      real  f3(D3_len)
      double precision  d3(D3_len)
      real  f11(D1_len, D1_len)
      double precision  d12(D2_len, D1_len)
      integer  s21(D1_len, D2_len)
      integer  i22(D2_len, D2_len)
      real  f23(D3_len, D2_len)
      integer  b32(D2_len, D3_len)
      integer  s33(D3_len, D3_len)
      integer  b112(D2_len, D1_len, D1_len)
      integer  s113(D3_len, D1_len, D1_len)
      real  f121(D1_len, D2_len, D1_len)
      double precision  d122(D2_len, D2_len, D1_len)
      integer  s131(D1_len, D3_len, D1_len)
      integer  i132(D2_len, D3_len, D1_len)
      real  f133(D3_len, D3_len, D1_len)
      real  f211(D1_len, D1_len, D2_len)
      double precision  d212(D2_len, D1_len, D2_len)
      integer  s221(D1_len, D2_len, D2_len)
      integer  i222(D2_len, D2_len, D2_len)
      real  f223(D3_len, D2_len, D2_len)
      integer  b232(D2_len, D3_len, D2_len)
      integer  s233(D3_len, D3_len, D2_len)
      integer  s311(D1_len, D1_len, D3_len)
      integer  i312(D2_len, D1_len, D3_len)
      real  f313(D3_len, D1_len, D3_len)
      double precision  var_dash_name_dash_dashes
      double precision  var_dot_name_dot_dots
* attribute vectors
      integer  int1val(4)
      integer  int2val(3)
      integer  intval(3)
      real  realval(3)
      double precision  doubleval(3)
* enter define mode
      iret = nf_create('ftest0.nc', NF_CLOBBER, ncid)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'Dr', NF_UNLIMITED, Dr_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'D1', 1, D1_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'D2', 2, D2_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'D3', 3, D3_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'dim-name-dashes', 4, dim_dash_name_dash_d
     1ashes_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'dim.name.dots', 5, dim_dot_name_dot_dots_
     1dim)
      call check_err(iret)
* define variables
      iret = nf_def_var(ncid, 'c', NF_CHAR, c_rank, 0, c_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'b', NF_INT1, b_rank, 0, b_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 's', NF_INT2, s_rank, 0, s_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'i', NF_INT, i_rank, 0, i_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'f', NF_REAL, f_rank, 0, f_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'd', NF_DOUBLE, d_rank, 0, d_id)
      call check_err(iret)
      cr_dims(1) = Dr_dim
      iret = nf_def_var(ncid, 'cr', NF_CHAR, cr_rank, cr_dims, cr_id)
      call check_err(iret)
      br_dims(1) = Dr_dim
      iret = nf_def_var(ncid, 'br', NF_INT1, br_rank, br_dims, br_id)
      call check_err(iret)
      sr_dims(1) = Dr_dim
      iret = nf_def_var(ncid, 'sr', NF_INT2, sr_rank, sr_dims, sr_id)
      call check_err(iret)
      ir_dims(1) = Dr_dim
      iret = nf_def_var(ncid, 'ir', NF_INT, ir_rank, ir_dims, ir_id)
      call check_err(iret)
      fr_dims(1) = Dr_dim
      iret = nf_def_var(ncid, 'fr', NF_REAL, fr_rank, fr_dims, fr_id)
      call check_err(iret)
      dr_dims(1) = Dr_dim
      iret = nf_def_var(ncid, 'dr', NF_DOUBLE, dr_rank, dr_dims, dr_id)
      call check_err(iret)
      c1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'c1', NF_CHAR, c1_rank, c1_dims, c1_id)
      call check_err(iret)
      b1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'b1', NF_INT1, b1_rank, b1_dims, b1_id)
      call check_err(iret)
      s1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 's1', NF_INT2, s1_rank, s1_dims, s1_id)
      call check_err(iret)
      i1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'i1', NF_INT, i1_rank, i1_dims, i1_id)
      call check_err(iret)
      f1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'f1', NF_REAL, f1_rank, f1_dims, f1_id)
      call check_err(iret)
      d1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'd1', NF_DOUBLE, d1_rank, d1_dims, d1_id)
      call check_err(iret)
      c2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'c2', NF_CHAR, c2_rank, c2_dims, c2_id)
      call check_err(iret)
      b2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'b2', NF_INT1, b2_rank, b2_dims, b2_id)
      call check_err(iret)
      s2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 's2', NF_INT2, s2_rank, s2_dims, s2_id)
      call check_err(iret)
      i2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'i2', NF_INT, i2_rank, i2_dims, i2_id)
      call check_err(iret)
      f2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'f2', NF_REAL, f2_rank, f2_dims, f2_id)
      call check_err(iret)
      d2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'd2', NF_DOUBLE, d2_rank, d2_dims, d2_id)
      call check_err(iret)
      c3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'c3', NF_CHAR, c3_rank, c3_dims, c3_id)
      call check_err(iret)
      b3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'b3', NF_INT1, b3_rank, b3_dims, b3_id)
      call check_err(iret)
      s3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 's3', NF_INT2, s3_rank, s3_dims, s3_id)
      call check_err(iret)
      i3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'i3', NF_INT, i3_rank, i3_dims, i3_id)
      call check_err(iret)
      f3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'f3', NF_REAL, f3_rank, f3_dims, f3_id)
      call check_err(iret)
      d3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'd3', NF_DOUBLE, d3_rank, d3_dims, d3_id)
      call check_err(iret)
      cr1_dims(2) = Dr_dim
      cr1_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'cr1', NF_CHAR, cr1_rank, cr1_dims, cr1_id
     1)
      call check_err(iret)
      br2_dims(2) = Dr_dim
      br2_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'br2', NF_INT1, br2_rank, br2_dims, br2_id
     1)
      call check_err(iret)
      sr3_dims(2) = Dr_dim
      sr3_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'sr3', NF_INT2, sr3_rank, sr3_dims, sr3_id
     1)
      call check_err(iret)
      f11_dims(2) = D1_dim
      f11_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'f11', NF_REAL, f11_rank, f11_dims, f11_id
     1)
      call check_err(iret)
      d12_dims(2) = D1_dim
      d12_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'd12', NF_DOUBLE, d12_rank, d12_dims, d12_
     1id)
      call check_err(iret)
      c13_dims(2) = D1_dim
      c13_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'c13', NF_CHAR, c13_rank, c13_dims, c13_id
     1)
      call check_err(iret)
      s21_dims(2) = D2_dim
      s21_dims(1) = D1_dim
      iret = nf_def_var(ncid, 's21', NF_INT2, s21_rank, s21_dims, s21_id
     1)
      call check_err(iret)
      i22_dims(2) = D2_dim
      i22_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'i22', NF_INT, i22_rank, i22_dims, i22_id)
      call check_err(iret)
      f23_dims(2) = D2_dim
      f23_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'f23', NF_REAL, f23_rank, f23_dims, f23_id
     1)
      call check_err(iret)
      c31_dims(2) = D3_dim
      c31_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'c31', NF_CHAR, c31_rank, c31_dims, c31_id
     1)
      call check_err(iret)
      b32_dims(2) = D3_dim
      b32_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'b32', NF_INT1, b32_rank, b32_dims, b32_id
     1)
      call check_err(iret)
      s33_dims(2) = D3_dim
      s33_dims(1) = D3_dim
      iret = nf_def_var(ncid, 's33', NF_INT2, s33_rank, s33_dims, s33_id
     1)
      call check_err(iret)
      sr11_dims(3) = Dr_dim
      sr11_dims(2) = D1_dim
      sr11_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'sr11', NF_INT2, sr11_rank, sr11_dims, sr1
     11_id)
      call check_err(iret)
      ir12_dims(3) = Dr_dim
      ir12_dims(2) = D1_dim
      ir12_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'ir12', NF_INT, ir12_rank, ir12_dims, ir12
     1_id)
      call check_err(iret)
      fr13_dims(3) = Dr_dim
      fr13_dims(2) = D1_dim
      fr13_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'fr13', NF_REAL, fr13_rank, fr13_dims, fr1
     13_id)
      call check_err(iret)
      cr21_dims(3) = Dr_dim
      cr21_dims(2) = D2_dim
      cr21_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'cr21', NF_CHAR, cr21_rank, cr21_dims, cr2
     11_id)
      call check_err(iret)
      br22_dims(3) = Dr_dim
      br22_dims(2) = D2_dim
      br22_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'br22', NF_INT1, br22_rank, br22_dims, br2
     12_id)
      call check_err(iret)
      sr23_dims(3) = Dr_dim
      sr23_dims(2) = D2_dim
      sr23_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'sr23', NF_INT2, sr23_rank, sr23_dims, sr2
     13_id)
      call check_err(iret)
      fr31_dims(3) = Dr_dim
      fr31_dims(2) = D3_dim
      fr31_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'fr31', NF_REAL, fr31_rank, fr31_dims, fr3
     11_id)
      call check_err(iret)
      dr32_dims(3) = Dr_dim
      dr32_dims(2) = D3_dim
      dr32_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'dr32', NF_DOUBLE, dr32_rank, dr32_dims, d
     1r32_id)
      call check_err(iret)
      cr33_dims(3) = Dr_dim
      cr33_dims(2) = D3_dim
      cr33_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'cr33', NF_CHAR, cr33_rank, cr33_dims, cr3
     13_id)
      call check_err(iret)
      c111_dims(3) = D1_dim
      c111_dims(2) = D1_dim
      c111_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'c111', NF_CHAR, c111_rank, c111_dims, c11
     11_id)
      call check_err(iret)
      b112_dims(3) = D1_dim
      b112_dims(2) = D1_dim
      b112_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'b112', NF_INT1, b112_rank, b112_dims, b11
     12_id)
      call check_err(iret)
      s113_dims(3) = D1_dim
      s113_dims(2) = D1_dim
      s113_dims(1) = D3_dim
      iret = nf_def_var(ncid, 's113', NF_INT2, s113_rank, s113_dims, s11
     13_id)
      call check_err(iret)
      f121_dims(3) = D1_dim
      f121_dims(2) = D2_dim
      f121_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'f121', NF_REAL, f121_rank, f121_dims, f12
     11_id)
      call check_err(iret)
      d122_dims(3) = D1_dim
      d122_dims(2) = D2_dim
      d122_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'd122', NF_DOUBLE, d122_rank, d122_dims, d
     1122_id)
      call check_err(iret)
      c123_dims(3) = D1_dim
      c123_dims(2) = D2_dim
      c123_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'c123', NF_CHAR, c123_rank, c123_dims, c12
     13_id)
      call check_err(iret)
      s131_dims(3) = D1_dim
      s131_dims(2) = D3_dim
      s131_dims(1) = D1_dim
      iret = nf_def_var(ncid, 's131', NF_INT2, s131_rank, s131_dims, s13
     11_id)
      call check_err(iret)
      i132_dims(3) = D1_dim
      i132_dims(2) = D3_dim
      i132_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'i132', NF_INT, i132_rank, i132_dims, i132
     1_id)
      call check_err(iret)
      f133_dims(3) = D1_dim
      f133_dims(2) = D3_dim
      f133_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'f133', NF_REAL, f133_rank, f133_dims, f13
     13_id)
      call check_err(iret)
      f211_dims(3) = D2_dim
      f211_dims(2) = D1_dim
      f211_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'f211', NF_REAL, f211_rank, f211_dims, f21
     11_id)
      call check_err(iret)
      d212_dims(3) = D2_dim
      d212_dims(2) = D1_dim
      d212_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'd212', NF_DOUBLE, d212_rank, d212_dims, d
     1212_id)
      call check_err(iret)
      c213_dims(3) = D2_dim
      c213_dims(2) = D1_dim
      c213_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'c213', NF_CHAR, c213_rank, c213_dims, c21
     13_id)
      call check_err(iret)
      s221_dims(3) = D2_dim
      s221_dims(2) = D2_dim
      s221_dims(1) = D1_dim
      iret = nf_def_var(ncid, 's221', NF_INT2, s221_rank, s221_dims, s22
     11_id)
      call check_err(iret)
      i222_dims(3) = D2_dim
      i222_dims(2) = D2_dim
      i222_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'i222', NF_INT, i222_rank, i222_dims, i222
     1_id)
      call check_err(iret)
      f223_dims(3) = D2_dim
      f223_dims(2) = D2_dim
      f223_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'f223', NF_REAL, f223_rank, f223_dims, f22
     13_id)
      call check_err(iret)
      c231_dims(3) = D2_dim
      c231_dims(2) = D3_dim
      c231_dims(1) = D1_dim
      iret = nf_def_var(ncid, 'c231', NF_CHAR, c231_rank, c231_dims, c23
     11_id)
      call check_err(iret)
      b232_dims(3) = D2_dim
      b232_dims(2) = D3_dim
      b232_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'b232', NF_INT1, b232_rank, b232_dims, b23
     12_id)
      call check_err(iret)
      s233_dims(3) = D2_dim
      s233_dims(2) = D3_dim
      s233_dims(1) = D3_dim
      iret = nf_def_var(ncid, 's233', NF_INT2, s233_rank, s233_dims, s23
     13_id)
      call check_err(iret)
      s311_dims(3) = D3_dim
      s311_dims(2) = D1_dim
      s311_dims(1) = D1_dim
      iret = nf_def_var(ncid, 's311', NF_INT2, s311_rank, s311_dims, s31
     11_id)
      call check_err(iret)
      i312_dims(3) = D3_dim
      i312_dims(2) = D1_dim
      i312_dims(1) = D2_dim
      iret = nf_def_var(ncid, 'i312', NF_INT, i312_rank, i312_dims, i312
     1_id)
      call check_err(iret)
      f313_dims(3) = D3_dim
      f313_dims(2) = D1_dim
      f313_dims(1) = D3_dim
      iret = nf_def_var(ncid, 'f313', NF_REAL, f313_rank, f313_dims, f31
     13_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'var-name-dashes', NF_DOUBLE, var_dash_nam
     1e_dash_dashes_rank, 0, var_dash_name_dash_dashes_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'var.name.dots', NF_DOUBLE, var_dot_name_d
     1ot_dots_rank, 0, var_dot_name_dot_dots_id)
      call check_err(iret)
* assign attributes
      intval(1) = 4
      iret = nf_put_att_int(ncid, c_id, 'att-name-dashes', NF_INT, 1, in
     1tval)
      call check_err(iret)
      intval(1) = 5
      iret = nf_put_att_int(ncid, c_id, 'att.name.dots', NF_INT, 1, intv
     1al)
      call check_err(iret)
      iret = nf_put_att_text(ncid, b_id, 'c', 1, char(0))
      call check_err(iret)
      int1val(1) = 0
      int1val(2) = 127
      int1val(3) = -128
      int1val(4) = -1
      iret = nf_put_att_int(ncid, s_id, 'b', NF_INT1, 4, int1val)
      call check_err(iret)
      int2val(1) = -32768
      int2val(2) = 0
      int2val(3) = 32767
      iret = nf_put_att_int(ncid, s_id, 's', NF_INT2, 3, int2val)
      call check_err(iret)
      intval(1) = -2147483647
      intval(2) = 0
      intval(3) = 2147483647
      iret = nf_put_att_int(ncid, i_id, 'i', NF_INT, 3, intval)
      call check_err(iret)
      realval(1) = -9.9999996e+35
      realval(2) = 0
      realval(3) = 9.9999996e+35
      iret = nf_put_att_real(ncid, i_id, 'f', NF_REAL, 3, realval)
      call check_err(iret)
      doubleval(1) = -1d+308
      doubleval(2) = 0
      doubleval(3) = 1d+308
      iret = nf_put_att_double(ncid, i_id, 'd', NF_DOUBLE, 3, doubleval)
      call check_err(iret)
      iret = nf_put_att_text(ncid, f_id, 'c', 1, 'x')
      call check_err(iret)
      iret = nf_put_att_text(ncid, d_id, 'c', 8, 'abcd'//char(9)//'Z$&')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'Gc', 1, char(0))
      call check_err(iret)
      int1val(1) = -128
      int1val(2) = 127
      iret = nf_put_att_int(ncid, NCGLOBAL, 'Gb', NF_INT1, 2, int1val)
      call check_err(iret)
      int2val(1) = -32768
      int2val(2) = 0
      int2val(3) = 32767
      iret = nf_put_att_int(ncid, NCGLOBAL, 'Gs', NF_INT2, 3, int2val)
      call check_err(iret)
      intval(1) = -2147483647
      intval(2) = 0
      intval(3) = 2147483647
      iret = nf_put_att_int(ncid, NCGLOBAL, 'Gi', NF_INT, 3, intval)
      call check_err(iret)
      realval(1) = -9.9999996e+35
      realval(2) = 0
      realval(3) = 9.9999996e+35
      iret = nf_put_att_real(ncid, NCGLOBAL, 'Gf', NF_REAL, 3, realval)
      call check_err(iret)
      doubleval(1) = -1d+308
      doubleval(2) = 0
      doubleval(3) = 1d+308
      iret = nf_put_att_double(ncid, NCGLOBAL, 'Gd', NF_DOUBLE, 3, doubl
     1eval)
      call check_err(iret)
      intval(1) = -1
      iret = nf_put_att_int(ncid, NCGLOBAL, 'Gatt-name-dashes', NF_INT, 
     11, intval)
      call check_err(iret)
      intval(1) = -2
      iret = nf_put_att_int(ncid, NCGLOBAL, 'Gatt.name.dots', NF_INT, 1,
     1 intval)
      call check_err(iret)
* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)
* store c
      iret = nf_put_var_text(ncid, c_id, '2')
      call check_err(iret)
* store b
      data b /-2/
      iret = nf_put_var_int(ncid, b_id, b)
      call check_err(iret)
* store s
      data s /-5/
      iret = nf_put_var_int(ncid, s_id, s)
      call check_err(iret)
* store i
      data i /-20/
      iret = nf_put_var_int(ncid, i_id, i)
      call check_err(iret)
* store f
      data f /-9/
      iret = nf_put_var_real(ncid, f_id, f)
      call check_err(iret)
* store d
      data d /-10./
      iret = nf_put_var_double(ncid, d_id, d)
      call check_err(iret)
* store c1
      iret = nf_put_var_text(ncid, c1_id, char(0))
      call check_err(iret)
* store b1
      data b1 /-128/
      iret = nf_put_var_int(ncid, b1_id, b1)
      call check_err(iret)
* store s1
      data s1 /-32768/
      iret = nf_put_var_int(ncid, s1_id, s1)
      call check_err(iret)
* store i1
      data i1 /-2147483646/
      iret = nf_put_var_int(ncid, i1_id, i1)
      call check_err(iret)
* store f1
      data f1 /-9.9999996e+35/
      iret = nf_put_var_real(ncid, f1_id, f1)
      call check_err(iret)
* store d1
      data d1 /-1.d+308/
      iret = nf_put_var_double(ncid, d1_id, d1)
      call check_err(iret)
* store c2
      iret = nf_put_var_text(ncid, c2_id, 'ab')
      call check_err(iret)
* store b2
      data b2 /-128, 127/
      iret = nf_put_var_int(ncid, b2_id, b2)
      call check_err(iret)
* store s2
      data s2 /-32768, 32767/
      iret = nf_put_var_int(ncid, s2_id, s2)
      call check_err(iret)
* store i2
      data i2 /-2147483646, 2147483647/
      iret = nf_put_var_int(ncid, i2_id, i2)
      call check_err(iret)
* store f2
      data f2 /-9.9999996e+35, 9.9999996e+35/
      iret = nf_put_var_real(ncid, f2_id, f2)
      call check_err(iret)
* store d2
      data d2 /-1.d+308, 1.d+308/
      iret = nf_put_var_double(ncid, d2_id, d2)
      call check_err(iret)
* store c3
      iret = nf_put_var_text(ncid, c3_id, char(1)//char(192)//'.')
      call check_err(iret)
* store b3
      data b3 /-128, 127, -1/
      iret = nf_put_var_int(ncid, b3_id, b3)
      call check_err(iret)
* store s3
      data s3 /-32768, 0, 32767/
      iret = nf_put_var_int(ncid, s3_id, s3)
      call check_err(iret)
* store i3
      data i3 /-2147483646, 0, 2147483647/
      iret = nf_put_var_int(ncid, i3_id, i3)
      call check_err(iret)
* store f3
      data f3 /-9.9999996e+35, 0, 9.9999996e+35/
      iret = nf_put_var_real(ncid, f3_id, f3)
      call check_err(iret)
* store d3
      data d3 /-1.d+308, 0., 1.d+308/
      iret = nf_put_var_double(ncid, d3_id, d3)
      call check_err(iret)
* store f11
      data f11 /-2187/
      iret = nf_put_var_real(ncid, f11_id, f11)
      call check_err(iret)
* store d12
      data d12 /-3000., -3010./
      iret = nf_put_var_double(ncid, d12_id, d12)
      call check_err(iret)
* store c13
      iret = nf_put_var_text(ncid, c13_id, char(9)//'b'//char(127))
      call check_err(iret)
* store s21
      data s21 /-375, -350/
      iret = nf_put_var_int(ncid, s21_id, s21)
      call check_err(iret)
* store i22
      data i22 /-24000, -24020, -23600, -23620/
      iret = nf_put_var_int(ncid, i22_id, i22)
      call check_err(iret)
* store f23
      data f23 /-2187, -2196, -2205, -2106, -2115, -2124/
      iret = nf_put_var_real(ncid, f23_id, f23)
      call check_err(iret)
* store c31
      iret = nf_put_var_text(ncid, c31_id, '+- ')
      call check_err(iret)
* store b32
      data b32 /-24, -26, -20, -22, -16, -18/
      iret = nf_put_var_int(ncid, b32_id, b32)
      call check_err(iret)
* store s33
      data s33 /-375, -380, -385, -350, -355, -360, -325, -330, -335/
      iret = nf_put_var_int(ncid, s33_id, s33)
      call check_err(iret)
* store c111
      iret = nf_put_var_text(ncid, c111_id, '@')
      call check_err(iret)
* store b112
      data b112 /64, 62/
      iret = nf_put_var_int(ncid, b112_id, b112)
      call check_err(iret)
* store s113
      data s113 /2500, 2495, 2490/
      iret = nf_put_var_int(ncid, s113_id, s113)
      call check_err(iret)
* store f121
      data f121 /26244, 26325/
      iret = nf_put_var_real(ncid, f121_id, f121)
      call check_err(iret)
* store d122
      data d122 /40000., 39990., 40100., 40090./
      iret = nf_put_var_double(ncid, d122_id, d122)
      call check_err(iret)
* store c123
      iret = nf_put_var_text(ncid, c123_id, 'one2'//char(0)//char(0))
      call check_err(iret)
* store s131
      data s131 /2500, 2525, 2550/
      iret = nf_put_var_int(ncid, s131_id, s131)
      call check_err(iret)
* store i132
      data i132 /640000, 639980, 640400, 640380, 640800, 640780/
      iret = nf_put_var_int(ncid, i132_id, i132)
      call check_err(iret)
* store f133
      data f133 /26244, 26235, 26226, 26325, 26316, 26307, 26406, 26397,
     1 26388/
      iret = nf_put_var_real(ncid, f133_id, f133)
      call check_err(iret)
* store f211
      data f211 /26244, 25515/
      iret = nf_put_var_real(ncid, f211_id, f211)
      call check_err(iret)
* store d212
      data d212 /40000., 39990., 39000., 38990./
      iret = nf_put_var_double(ncid, d212_id, d212)
      call check_err(iret)
* store s221
      data s221 /2500, 2525, 2375, 2400/
      iret = nf_put_var_int(ncid, s221_id, s221)
      call check_err(iret)
* store i222
      data i222 /640000, 639980, 640400, 640380, 632000, 631980, 632400,
     1 632380/
      iret = nf_put_var_int(ncid, i222_id, i222)
      call check_err(iret)
* store f223
      data f223 /26244, 26235, 26226, 26325, 26316, 26307, 25515, 25506,
     1 25497, 25596, 25587, 25578/
      iret = nf_put_var_real(ncid, f223_id, f223)
      call check_err(iret)
* store c231
      iret = nf_put_var_text(ncid, c231_id, '@DHHLP')
      call check_err(iret)
* store b232
      data b232 /64, 62, 68, 66, 72, 70, 56, 54, 60, 58, 64, 62/
      iret = nf_put_var_int(ncid, b232_id, b232)
      call check_err(iret)
* store s233
      data s233 /2500, 2495, 2490, 2525, 2520, 2515, 2550, 2545, 2540, 2
     1375, 2370, 2365, 2400, 2395, 2390, 2425, 2420, 2415/
      iret = nf_put_var_int(ncid, s233_id, s233)
      call check_err(iret)
* store s311
      data s311 /2500, 2375, 2250/
      iret = nf_put_var_int(ncid, s311_id, s311)
      call check_err(iret)
* store i312
      data i312 /640000, 639980, 632000, 631980, 624000, 623980/
      iret = nf_put_var_int(ncid, i312_id, i312)
      call check_err(iret)
* store f313
      data f313 /26244, 26235, 26226, 25515, 25506, 25497, 24786, 24777,
     1 24768/
      iret = nf_put_var_real(ncid, f313_id, f313)
      call check_err(iret)
* store var-name-dashes
      data var_dash_name_dash_dashes /-1./
      iret = nf_put_var_double(ncid, var_dash_name_dash_dashes_id, var_d
     1ash_name_dash_dashes)
      call check_err(iret)
* store var.name.dots
      data var_dot_name_dot_dots /-2./
      iret = nf_put_var_double(ncid, var_dot_name_dot_dots_id, var_dot_n
     1ame_dot_dots)
      call check_err(iret)
       
* Write record variables
      call writerecs(ncid,cr_id,br_id,sr_id,ir_id,fr_id,dr_id,cr1_id,br2
     1_id,sr3_id,sr11_id,ir12_id,fr13_id,cr21_id,br22_id,sr23_id,fr31_id
     2,dr32_id,cr33_id)
       
      iret = nf_close(ncid)
      call check_err(iret)
      end
       
      subroutine writerecs(ncid,cr_id,br_id,sr_id,ir_id,fr_id,dr_id,cr1_
     1id,br2_id,sr3_id,sr11_id,ir12_id,fr13_id,cr21_id,br22_id,sr23_id,f
     2r31_id,dr32_id,cr33_id)
       
* netCDF id
      integer  ncid
* variable ids
      integer  cr_id
      integer  br_id
      integer  sr_id
      integer  ir_id
      integer  fr_id
      integer  dr_id
      integer  cr1_id
      integer  br2_id
      integer  sr3_id
      integer  sr11_id
      integer  ir12_id
      integer  fr13_id
      integer  cr21_id
      integer  br22_id
      integer  sr23_id
      integer  fr31_id
      integer  dr32_id
      integer  cr33_id
       
      include 'netcdf.inc'
* error status return
      integer  iret
       
* netCDF dimension sizes for dimensions used with record variables
      integer  D1_len
      parameter (D1_len = 1)
      integer  D2_len
      parameter (D2_len = 2)
      integer  D3_len
      parameter (D3_len = 3)
       
* rank (number of dimensions) for each variable
      integer  cr_rank
      integer  br_rank
      integer  sr_rank
      integer  ir_rank
      integer  fr_rank
      integer  dr_rank
      integer  cr1_rank
      integer  br2_rank
      integer  sr3_rank
      integer  sr11_rank
      integer  ir12_rank
      integer  fr13_rank
      integer  cr21_rank
      integer  br22_rank
      integer  sr23_rank
      integer  fr31_rank
      integer  dr32_rank
      integer  cr33_rank
      parameter (cr_rank = 1)
      parameter (br_rank = 1)
      parameter (sr_rank = 1)
      parameter (ir_rank = 1)
      parameter (fr_rank = 1)
      parameter (dr_rank = 1)
      parameter (cr1_rank = 2)
      parameter (br2_rank = 2)
      parameter (sr3_rank = 2)
      parameter (sr11_rank = 3)
      parameter (ir12_rank = 3)
      parameter (fr13_rank = 3)
      parameter (cr21_rank = 3)
      parameter (br22_rank = 3)
      parameter (sr23_rank = 3)
      parameter (fr31_rank = 3)
      parameter (dr32_rank = 3)
      parameter (cr33_rank = 3)
* starts and counts for array sections of record variables
      integer  cr_start(cr_rank), cr_count(cr_rank)
      integer  br_start(br_rank), br_count(br_rank)
      integer  sr_start(sr_rank), sr_count(sr_rank)
      integer  ir_start(ir_rank), ir_count(ir_rank)
      integer  fr_start(fr_rank), fr_count(fr_rank)
      integer  dr_start(dr_rank), dr_count(dr_rank)
      integer  cr1_start(cr1_rank), cr1_count(cr1_rank)
      integer  br2_start(br2_rank), br2_count(br2_rank)
      integer  sr3_start(sr3_rank), sr3_count(sr3_rank)
      integer  sr11_start(sr11_rank), sr11_count(sr11_rank)
      integer  ir12_start(ir12_rank), ir12_count(ir12_rank)
      integer  fr13_start(fr13_rank), fr13_count(fr13_rank)
      integer  cr21_start(cr21_rank), cr21_count(cr21_rank)
      integer  br22_start(br22_rank), br22_count(br22_rank)
      integer  sr23_start(sr23_rank), sr23_count(sr23_rank)
      integer  fr31_start(fr31_rank), fr31_count(fr31_rank)
      integer  dr32_start(dr32_rank), dr32_count(dr32_rank)
      integer  cr33_start(cr33_rank), cr33_count(cr33_rank)
       
* data variables
       
      integer  cr_nr
      parameter (cr_nr = 2)
       
      integer  br_nr
      parameter (br_nr = 2)
      integer  br(br_nr)
       
      integer  sr_nr
      parameter (sr_nr = 2)
      integer  sr(sr_nr)
       
      integer  ir_nr
      parameter (ir_nr = 2)
      integer  ir(ir_nr)
       
      integer  fr_nr
      parameter (fr_nr = 2)
      real  fr(fr_nr)
       
      integer  dr_nr
      parameter (dr_nr = 2)
      double precision  dr(dr_nr)
       
      integer  cr1_nr
      parameter (cr1_nr = 2)
       
      integer  br2_nr
      parameter (br2_nr = 2)
      integer  br2(D2_len, br2_nr)
       
      integer  sr3_nr
      parameter (sr3_nr = 2)
      integer  sr3(D3_len, sr3_nr)
       
      integer  sr11_nr
      parameter (sr11_nr = 2)
      integer  sr11(D1_len, D1_len, sr11_nr)
       
      integer  ir12_nr
      parameter (ir12_nr = 2)
      integer  ir12(D2_len, D1_len, ir12_nr)
       
      integer  fr13_nr
      parameter (fr13_nr = 2)
      real  fr13(D3_len, D1_len, fr13_nr)
       
      integer  cr21_nr
      parameter (cr21_nr = 2)
       
      integer  br22_nr
      parameter (br22_nr = 2)
      integer  br22(D2_len, D2_len, br22_nr)
       
      integer  sr23_nr
      parameter (sr23_nr = 2)
      integer  sr23(D3_len, D2_len, sr23_nr)
       
      integer  fr31_nr
      parameter (fr31_nr = 2)
      real  fr31(D1_len, D3_len, fr31_nr)
       
      integer  dr32_nr
      parameter (dr32_nr = 2)
      double precision  dr32(D2_len, D3_len, dr32_nr)
       
      integer  cr33_nr
      parameter (cr33_nr = 2)
       
      data br /-128, 127/
      data sr /-32768, 32767/
      data ir /-2147483646, 2147483647/
      data fr /-9.9999996e+35, 9.9999996e+35/
      data dr /-1.d+308, 1.d+308/
      data br2 /-24, -26, -20, -22/
      data sr3 /-375, -380, -385, -350, -355, -360/
      data sr11 /2500, 2375/
      data ir12 /640000, 639980, 632000, 631980/
      data fr13 /26244, 26235, 26226, 25515, 25506, 25497/
      data br22 /64, 62, 68, 66, 56, 54, 60, 58/
      data sr23 /2500, 2495, 2490, 2525, 2520, 2515, 2375, 2370, 2365, 2
     1400, 2395, 2390/
      data fr31 /26244, 26325, 26406, 25515, 25596, 25677/
      data dr32 /40000., 39990., 40100., 40090., 40200., 40190., 39000.,
     1 38990., 39100., 39090., 39200., 39190./
       
* store cr
      cr_start(1) = 1
      cr_count(1) = cr_nr
      iret = nf_put_vara_text(ncid, cr_id, cr_start, cr_count, 'ab')
      call check_err(iret)
* store br
      br_start(1) = 1
      br_count(1) = br_nr
      iret = nf_put_vara_int(ncid, br_id, br_start, br_count, br)
      call check_err(iret)
* store sr
      sr_start(1) = 1
      sr_count(1) = sr_nr
      iret = nf_put_vara_int(ncid, sr_id, sr_start, sr_count, sr)
      call check_err(iret)
* store ir
      ir_start(1) = 1
      ir_count(1) = ir_nr
      iret = nf_put_vara_int(ncid, ir_id, ir_start, ir_count, ir)
      call check_err(iret)
* store fr
      fr_start(1) = 1
      fr_count(1) = fr_nr
      iret = nf_put_vara_real(ncid, fr_id, fr_start, fr_count, fr)
      call check_err(iret)
* store dr
      dr_start(1) = 1
      dr_count(1) = dr_nr
      iret = nf_put_vara_double(ncid, dr_id, dr_start, dr_count, dr)
      call check_err(iret)
* store cr1
      cr1_start(1) = 1
      cr1_start(2) = 1
      cr1_count(1) = D1_len
      cr1_count(2) = cr1_nr
      iret = nf_put_vara_text(ncid, cr1_id, cr1_start, cr1_count, 'xy')
      call check_err(iret)
* store br2
      br2_start(1) = 1
      br2_start(2) = 1
      br2_count(1) = D2_len
      br2_count(2) = br2_nr
      iret = nf_put_vara_int(ncid, br2_id, br2_start, br2_count, br2)
      call check_err(iret)
* store sr3
      sr3_start(1) = 1
      sr3_start(2) = 1
      sr3_count(1) = D3_len
      sr3_count(2) = sr3_nr
      iret = nf_put_vara_int(ncid, sr3_id, sr3_start, sr3_count, sr3)
      call check_err(iret)
* store sr11
      sr11_start(1) = 1
      sr11_start(2) = 1
      sr11_start(3) = 1
      sr11_count(1) = D1_len
      sr11_count(2) = D1_len
      sr11_count(3) = sr11_nr
      iret = nf_put_vara_int(ncid, sr11_id, sr11_start, sr11_count, sr11
     1)
      call check_err(iret)
* store ir12
      ir12_start(1) = 1
      ir12_start(2) = 1
      ir12_start(3) = 1
      ir12_count(1) = D2_len
      ir12_count(2) = D1_len
      ir12_count(3) = ir12_nr
      iret = nf_put_vara_int(ncid, ir12_id, ir12_start, ir12_count, ir12
     1)
      call check_err(iret)
* store fr13
      fr13_start(1) = 1
      fr13_start(2) = 1
      fr13_start(3) = 1
      fr13_count(1) = D3_len
      fr13_count(2) = D1_len
      fr13_count(3) = fr13_nr
      iret = nf_put_vara_real(ncid, fr13_id, fr13_start, fr13_count, fr1
     13)
      call check_err(iret)
* store cr21
      cr21_start(1) = 1
      cr21_start(2) = 1
      cr21_start(3) = 1
      cr21_count(1) = D1_len
      cr21_count(2) = D2_len
      cr21_count(3) = cr21_nr
      iret = nf_put_vara_text(ncid, cr21_id, cr21_start, cr21_count, '@D
     1HL')
      call check_err(iret)
* store br22
      br22_start(1) = 1
      br22_start(2) = 1
      br22_start(3) = 1
      br22_count(1) = D2_len
      br22_count(2) = D2_len
      br22_count(3) = br22_nr
      iret = nf_put_vara_int(ncid, br22_id, br22_start, br22_count, br22
     1)
      call check_err(iret)
* store sr23
      sr23_start(1) = 1
      sr23_start(2) = 1
      sr23_start(3) = 1
      sr23_count(1) = D3_len
      sr23_count(2) = D2_len
      sr23_count(3) = sr23_nr
      iret = nf_put_vara_int(ncid, sr23_id, sr23_start, sr23_count, sr23
     1)
      call check_err(iret)
* store fr31
      fr31_start(1) = 1
      fr31_start(2) = 1
      fr31_start(3) = 1
      fr31_count(1) = D1_len
      fr31_count(2) = D3_len
      fr31_count(3) = fr31_nr
      iret = nf_put_vara_real(ncid, fr31_id, fr31_start, fr31_count, fr3
     11)
      call check_err(iret)
* store dr32
      dr32_start(1) = 1
      dr32_start(2) = 1
      dr32_start(3) = 1
      dr32_count(1) = D2_len
      dr32_count(2) = D3_len
      dr32_count(3) = dr32_nr
      iret = nf_put_vara_double(ncid, dr32_id, dr32_start, dr32_count, d
     1r32)
      call check_err(iret)
* store cr33
      cr33_start(1) = 1
      cr33_start(2) = 1
      cr33_start(3) = 1
      cr33_count(1) = D3_len
      cr33_count(2) = D3_len
      cr33_count(3) = cr33_nr
      iret = nf_put_vara_text(ncid, cr33_id, cr33_start, cr33_count, '1'
     1//char(0)//char(0)//'two3'//char(0)//char(0)//'4'//char(0)//char(0
     2)//'5'//char(0)//char(0)//'six')
      call check_err(iret)
       
      end
       
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end
