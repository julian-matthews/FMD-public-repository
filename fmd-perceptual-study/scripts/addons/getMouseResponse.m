%Mouse response
function [x,y,buttons] = getMouseResponse()

[x,y,buttons] = GetMouse;
while any(buttons) % if already down, wait for release
    [x,y,buttons] = GetMouse;
end
while ~any(buttons) % wait for press
    [x,y,buttons] = GetMouse;
end
end

%Check inside image
function rect = checkRectInsideImage (image, rect);

if rect(1) < 0
    rect(1) = 0;
end

if rect(3) > size(image,2)
    rect(3) = size(image,2);
end

if rect(2) < 0
    rect(2) = 0;
end

if rect(4) > size(image, 1)
    rect(4) = size(image, 1);
end
end

