function mObject = updateFieldOfMobject(mObject,fieldName,fieldValue)

if isfield(mObject,fieldName)
    mObject = setfield(mObject,fieldName,fieldValue);
else
    disp([fieldName ' is not in a Mobject.']);
end;
    
    