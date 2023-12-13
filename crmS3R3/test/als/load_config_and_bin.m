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

maxtimesteps = 15*2^18;

% First, read the lua file
luafname = [fname '.lua'];
try
out = readLua(luafname, ...
   {... Integrator
    'integrator',...
    ... Integrator options
    'const_mass_matrix',...
    'diag_mass_matrix',...
    'banded_iteration_matrix',...
    'nr_subdiag',...
    'nr_superdiag',...
    'recalc_iteration_matrix',...
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
    'problem_name',...
    'kirchhoff',...
    'inextensible',...
    'use_SMA',...
    'n',...
    'fixed_x0',...
    'fixed_x0_position',...
    'fixed_xn',...
    'fixed_xn_position',...
    'fixed_p0',...
    'fixed_p0_orientation',...
    'fixed_pn',...
    'fixed_pn_orientation',...
    'L',...
    'ds',...
    'mI',...
    'm',...
    'CGam',...
    'CGamd',...
    'CK',...
    'CKd',...
    ... Output options
    'output_s_at',...
    'output_s',...
    'output_t_at',...
    't_output_at_multiples_of',...
    ... External forces and moments
    'external',...
    'variable_steps'});
   if strcmp(out.integrator,'gena') || isempty(out.integrator)
      intopts = readLua(luafname, ...
         {'rho_infinity',...
          'alpha_m',...
          'alpha_f',...
          'beta',...
          'gamma',...
          'perturb',...
          'perturb_s'});
   elseif strcmp(out.integrator,'BLieDF')
      intopts = readLua(luafname, {'k_bdf'});
   else
      intopts = struct();
   end

   % Merge integrator options with out
   f = fieldnames(intopts);
   for j = 1:length(f)
      out.(f{j}) = intopts.(f{j});
   end
catch ME
   warning(['Error reading Lua file, skipping' char(10) ME.message]);
   out = [];
   matched = false;
   return;
end

if isempty(out.steps)
   warning([luafname ' has no ''steps'', so it''s probably empty or corrupted. Skipping']);
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

disp(luafname); % DEBUG

%% Calculate sizes
if strcmp(out.integrator,'RATTLie') || strcmp(out.integrator,'SHAKELie')
   has_vd = 0;
else
   has_vd = 1;
end
if (out.output_s_at == 1)
   sizeq = 7*length(out.output_s);
   sizev = 6*length(out.output_s);
   sizel = 0;
   if out.kirchhoff == 1
      sizel = sizel + 2*length(out.output_s);
   end
%   if out.inextensible == 1
%      sizel = sizel + length(out.s_output_p);
%   end
   sizebin1 = 1 + sizeq + (1+has_vd)*sizev + sizel;
   if (out.stab2 == 1)
      sizebin1 = sizebin1 + sizel;
   end
else
   sizeq = (4+3)*(out.n+1);
   sizev = (3+3)*(out.n+1);
   if (out.fixed_x0 == 1)
      sizeq = sizeq - 3;
      sizev = sizev - 3;
   end
   if (out.fixed_xn == 1)
      sizeq = sizeq - 3;
      sizev = sizev - 3;
   end
   if (out.fixed_p0 == 1)
      sizeq = sizeq - 4;
      sizev = sizev - 3;
   end
   if (out.fixed_pn == 1)
      sizeq = sizeq - 4;
      sizev = sizev - 3;
   end
   sizel = 0;
   if out.kirchhoff == 1
      sizel = sizel + 2*out.n;
   end
%   if out.inextensible == 1
%      sizel = sizel + (out.n - 1);
%   end

   sizebin1 = 1 + sizeq + (1+has_vd)*sizev + 3*sizel;
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
   if isempty(out.variable_steps) || out.variable_steps ~= 1
      sizebin2 = out.steps + 1;
   else
      sizebin2 = 4*out.steps + 1
   end
end



% Test if binary file is intact
binfname = [fname '.bin'];
if isempty(out.variable_steps) || out.variable_steps ~= 1
   dr = dir(binfname);
   if dr.bytes ~= sizebin1 * sizebin2 * 8
      warning([binfname ' is not complete (or corrupted)']);
      out.rslt.finished = false;
   end
end


% Next, read the binary file
binfhandle = fopen(binfname);
bin = fread(binfhandle, [sizebin1, min(sizebin2,maxtimesteps)], 'real*8');
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
if has_vd
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
   if not(out.output_s_at == 1)
      out.rslt.Phi = bin(1+sizeq+(1+has_vd)*sizev+sizel+sizeeta + 1:...
                         1+sizeq+(1+has_vd)*sizev+sizel+sizeeta + sizel,:);
      out.rslt.Bv  = bin(1+sizeq+(1+has_vd)*sizev+2*sizel+sizeeta + 1:...
                         1+sizeq+(1+has_vd)*sizev+2*sizel+sizeeta + sizel,:);
   end
end

% Check second dimension of bin
if isempty(out.variable_steps) || out.variable_steps ~= 1
   if size(bin,2) ~= sizebin2
      out.rslt.finished = false;
      return
   elseif size(bin,2) == maxtimesteps
      warning(['Only did read up to ' num2str(maxtimesteps) ' time steps (maxtimesteps)']);
      out.rslt.finished = false;
   else
      % Everything succeeded
      out.rslt.finished = true;
   end
else
   if abs(out.rslt.t(end) - out.te) > 1e-5
      out.rslt.finished = false;
      warning([binfname ' is not complete (or corrupted)']);
      return
   else
      out.rslt.finished = true;
   end
end

% Load statistics from lua file
out.stats = readLua(luafname, ...
                    {'cpu_time',...
                     'newt_steps_max',...
                     'newt_steps_avg',...
                     'n_g_calls',...
                     'n_B_calls'});
