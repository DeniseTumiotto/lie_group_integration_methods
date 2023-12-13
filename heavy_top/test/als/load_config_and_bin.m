function [out, matched] = load_config_and_bin(fname, pattern, only_conf)
% Loads the lua file [fname '.lua'] and the binary file [fname .'bin'].
%
% If the second argument, a struct pattern, is given this function will
% check if the struct obtained from the lua file matches it. If not, the
% binary file will not be read and the second return argument will be
% set to false.
%
% If the third argument is given and set to 'only_conf', the binary will
% will not be loaded.

% First, read the lua file
luafname = [fname '.lua'];
out = readLua(luafname, ...
   {... Integrator
    'integrator',...
    ... Problemset
    'problemset',...
    ... Algorithmic parameters
    'alpha_m',...
    'alpha_f',...
    'beta',...
    'gamma',...
    ... Integrator options
    'const_mass_matrix',...
    'diag_mass_matrix',...
    'banded_iteration_matrix',...
    'nr_subdiag',...
    'nr_superdiag',...
    'recalc_iteration_matrix',...
    'perturb',...
    'perturb_s',...
    'use_num_K',...
    'use_num_D',...
    'no_K',...
    'no_D',...
    'rtol',...
    'atol',...
    'imax',...
    'stab2',...
    ... Integration interval and step size
    't0',...
    'te',...
    'steps',...
    ... Problem options
    'mass',...
    'gravity',...
    'inerJ',...
    'x0',...
    'p0',...
    'Om0',...
    'liegroup',...
    ... Output options
    'output_type',...
    'output_t_at',...
    't_output_at_multiples_of'});

if isempty(out.alpha_m)
   warning([luafname ' has no alpha_m, so it''s probably empty or corrupted. Skipping']);
   matched = false;
   return;
end

if nargin >= 2 && ~structmatch(out, pattern)
   matched = false;
   return;
else
   matched = true;
end

if nargin >= 3 && strcmp(only_conf,'only_conf')
   return;
end

% Calculate sizes
if strcmp(out.integrator,'RATTLie') || strcmp(out.integrator,'SHAKELie')
   has_vd = 0;
else
   has_vd = 1;
end
switch out.liegroup
   case 1 % SO(3)
      sizeq = 9;
      sizev = 3;
      sizel = -1;
   case 2 % S³
      sizeq = 4;
      sizev = 3;
      sizel = -1;
   case 3 % SO(3) x R³
      sizeq = 9+3;
      sizev = 3+3;
      sizel = 3;
   case 4 % S³ x R³
      sizeq = 4+3;
      sizev = 3+3;
      sizel = 3;
   case 5 % SE(3)
      sizeq = 9+3;
      sizev = 3+3;
      sizel = 3;
   case 6 % S³ |x R³
      sizeq = 4+3;
      sizev = 3+3;
      sizel = 3;
   case 7 % UDQ
      sizeq = 4+4;
      sizev = 3+3;
      sizel = 3;
end
sizebin1 = 1 + sizeq + (1 + has_vd)*sizev;
if sizel > 0
   sizebin1 = sizebin1 + 3*sizel;
   if (out.stab2 == 1)
      sizebin1 = sizebin1 + sizel;
   end
end

if (out.output_t_at == 1)
   if not(out.t0 == 0.0)
      error('output_t_at == 1, but not t0 == 0.0');
   end
   sizebin2 = floor(out.te/out.t_output_at_multiples_of) + 1;
else
   sizebin2 = out.steps + 1;
end

% Test if binary file is intact
binfname = [fname '.bin'];
dr = dir(binfname);
if dr.bytes ~= sizebin1 * sizebin2 * 8
   warning([binfname 'is not complete (or corrupted)']);
   out.rslt.finished = false;
end

% Next, read the binary file
binfhandle = fopen(binfname);
bin = fread(binfhandle, [sizebin1, sizebin2], 'real*8');
fclose(binfhandle);

% Check if the first dimension of bin agrees
if size(bin,1) ~= sizebin1
   out.rslt.finished = false;
   return;
end

% Add q, v, vd etc. to out
out.rslt.t  = bin(1,:);
out.rslt.q  = bin(1 + 1:...
                  1 + sizeq,:);
out.rslt.v  = bin(1+sizeq + 1:...
                  1+sizeq + sizev,:);
if has_vd ~= 0
out.rslt.vd = bin(1+sizeq+sizev + 1:...
                  1+sizeq+sizev + sizev,:);
end
if (sizel > 0)
   out.rslt.l = bin(1+sizeq+(1+has_vd)*sizev + 1:...
                    1+sizeq+(1+has_vd)*sizev + sizel,:);
   if (out.stab2 == 1)
      out.rslt.e = bin(1+sizeq+(1+has_vd)*sizev+sizel + 1:...
                       1+sizeq+(1+has_vd)*sizev+sizel + sizel,:);
      sizeeta = sizel;
   else
      sizeeta = 0;
   end
   out.rslt.Phi = bin(1+sizeq+(1+has_vd)*sizev+sizel+sizeeta + 1:...
                      1+sizeq+(1+has_vd)*sizev+sizel+sizeeta + sizel,:);
   out.rslt.Bv  = bin(1+sizeq+(1+has_vd)*sizev+2*sizel+sizeeta + 1:...
                      1+sizeq+(1+has_vd)*sizev+2*sizel+sizeeta + sizel,:);
end

% Check second dimension of bin
if size(bin,2) ~= sizebin2
   out.rslt.finished = false;
   return
else
   % Everything succeeded
   out.rslt.finished = true;
end

% Load statistics from lua file
out.stats = readLua(luafname, ...
                    {'cpu_time',...
                     'newt_steps_max',...
                     'newt_steps_avg',...
                     'n_g_calls',...
                     'n_B_calls'});
