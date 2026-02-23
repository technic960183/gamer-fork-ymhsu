import pandas as pd
import matplotlib.pyplot as plt

base_name = 'Xline_y0.000_z0.000_'
file_range = [0, 1, 2]
# file_range = list(range(0, 3))
# file_range = [0, 5]

plt.figure(figsize=(4, 6))

for n in file_range:
    file_path = f'{base_name}{n:06d}'
    with open(file_path, 'r') as f:
        lines = f.readlines()

    header = [h for h in lines[0].strip().split()[1:] if h]
    data_lines = lines[1:]
    data = [list(map(float, line.strip().split())) for line in data_lines]
    df = pd.DataFrame(data, columns=header)
    # df.set_index('i', inplace=True)

    plt.plot(df['x'], df['CR_E'], label=f't={n*0.01:.2f}')

plt.xlabel('x')
plt.ylabel('CR_E')
plt.ylim(0.9, 2.1)
plt.title('GAMER CR_E vs x')
plt.legend(title='File')
plt.grid(True)
plt.tight_layout()
plt.savefig('CR_E_vs_i.png')
