%'Odometer pattern'

possible_values = {'Z_Z', -3, 'V7', 'G', 'C', 119, 'Sr', 'Q', '2J'};
num_diff_vals = length(possible_values);
NumPlaces = 5;   %or as needed
num_possibilities = num_diff_vals^NumPlaces;
results = cell(num_possibilities, NumPlaces);
odo = zeros(1, NumPlaces);
for K = 1 : num_possibilities
  results(K,:) = possible_values( odo + 1 );    %note: adding 1 to each position
  for P = NumPlaces : -1 : 1
    odo(P) = odo(P) + 1;
    if odo(P) == num_diff_vals
      odo(P) = 0;
    else
      break;
    end
  end
end
