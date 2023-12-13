function out = structmatch(struct, pattern)
% Tests if a struct matches a pattern struct.
%
% The struct struct matches the struct pattern, if all fields of pattern
% exists in struct and have the same values.
%
% Note that pattern may only contain fields that are numeric, logical,
% strings, cell array of strings or structs, otherwise structmatch will
% throw an error.

fields = fieldnames(pattern);

for i=1:numel(fields)
   field = fields{i};
   val = pattern.(field);
   
   if isnumeric(val) || islogical(val)
      test = (struct.(field) ~= val);
      if any(test(:));
         out = false;
         return;
      end
   elseif ischar(val) || iscellstr(val)
      if ~strcmp(struct.(field), val)
         out = false;
         return;
      end
   elseif isstruct(val)
      if ~isfield(struct,field) || ~structmatch(struct.(field), val)
         out = false;
         return;
      end
   else
      error('structmatch:unsupportedClass',...
         ['Struct field ' field ' has unsupported class ' class(val)]);
   end
end

out = true;