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

only_nonhol = 0;

% First, read the lua file
luafname = [fname '.lua'];
try
out = readLua(luafname, ...
   {... Integrator options
    'const_mass_matrix',...
    'diag_mass_matrix',...
    'banded_iteration_matrix',...
    'nr_subdiag',...
    'nr_superdiag',...
    'recalc_iteration_matrix',...
    'rtol',...
    'atol',...
    'imax',...
    ... Integration interval and step size
    't0',...
    'te',...
    'steps',...
    ... Problem options
    'problem_name',...
    'mass',...
    'gravity',...
    'inerJ',...
    'x0',...
    'p0',...
    'Om0',...
    ... Output options
    'output_t_at',...
    't_output_at_multiples_of'});
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
has_vd = 0;
sizeq = 7;
sizev = 6;
if only_nonhol
   sizel = 0;
   sizelnh = 3;
else
   sizel = 1;
   sizelnh = 2;
end
%          t   q       v (vd)             l lm lnh Phi Bv Bnhv
sizebin1 = 1 + sizeq + (1+has_vd)*sizev + 4*sizel + 2*sizelnh;
%
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
   warning([binfname ' is not complete (or corrupted)']);
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
if has_vd               
   out.rslt.vd = bin(1+sizeq+sizev + 1:...
                     1+sizeq+sizev + sizev,:);
end                  
if (sizel > 0)
   out.rslt.l = bin(1+sizeq+(1+has_vd)*sizev + 1:...
                    1+sizeq+(1+has_vd)*sizev + sizel,:);
   out.rslt.lm = bin(1+sizeq+(1+has_vd)*sizev+sizel + 1:...
                     1+sizeq+(1+has_vd)*sizev+sizel + sizel,:);
   out.rslt.lnh = bin(1+sizeq+(1+has_vd)*sizev+2*sizel + 1:...
                      1+sizeq+(1+has_vd)*sizev+2*sizel + sizelnh,:);                  
   out.rslt.Phi = bin(1+sizeq+(1+has_vd)*sizev+2*sizel+sizelnh + 1:...
                      1+sizeq+(1+has_vd)*sizev+3*sizel+sizelnh,:);
   out.rslt.Bv  = bin(1+sizeq+(1+has_vd)*sizev+3*sizel+sizelnh + 1:...
                      1+sizeq+(1+has_vd)*sizev+4*sizel+sizelnh,:);
   out.rslt.Bnhv= bin(1+sizeq+(1+has_vd)*sizev+4*sizel+sizelnh + 1:...
                      1+sizeq+(1+has_vd)*sizev+4*sizel+2*sizelnh,:);                   
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