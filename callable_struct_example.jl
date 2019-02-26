mutable struct AddFunction
  a
end

function (f::AddFunction)(x)
  f.a + x
end

function AddFunction()
  AddFunction(5)
end

function AddFunction(x,y,z)
  AddFunction(x + y - z)
end
