function bool = isrealsymmatrix(M)

bool = 1;

sz = size(M);

if sz(1) ~= sz(2)
  bool = 0;
  return
end

if ~isnumeric(M) || ~isreal(M)
  bool = 0;
end



if ~all(all(abs(M - M') < 1e-12))
%  error('expecting a symmetric matrix in isrealsymmatrix')
  bool = 0;
end
  
