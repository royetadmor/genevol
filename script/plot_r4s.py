import re
import matplotlib.pyplot as plt

sites, rates = [], []
with open("r4s_data_9dictos_100.txt") as f:
    for line in f:
        m = re.match(r"Expected rate at site (\d+) is ([0-9.]+)", line.strip())
        if m:
            sites.append(int(m.group(1)))
            rates.append(float(m.group(2)))

plt.figure(figsize=(12, 4))
plt.plot(sites, rates, linewidth=0.8, color="steelblue")
plt.axhline(y=1, color="red", linewidth=1.0, linestyle="--", label="Mean (1.0)")
plt.legend()
plt.xlabel("Site")
plt.ylabel("Expected rate")
plt.title("Expected rate per site (15 species, 1000 families)")
plt.tight_layout()
plt.savefig("r4s_data_Eudicots.png", dpi=150)
print("Saved Eudicots.png")

s = 0
for r in rates:
    s += r
s /= len(rates)
print("avg is ", s)
