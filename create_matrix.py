import random

names = 'A B C D E F G H I K L M N O P Q R S T V X Y Z'.split()
count = 20

with open('tests/test5.txt', 'w') as f:
	f.write(' '.join(names[:count]))
	f.write('\n')

	for i in range(count):
		distances = [names[i]] + ['0'] * (i + 1) + [str(10*random.random()) for _ in range(count - i - 1)]
		f.write(' '.join(distances))
		f.write('\n')
