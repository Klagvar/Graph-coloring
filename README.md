# Graph-coloring
## Описание
Этот проект представляет собой реализацию шести приближенных алгоритмов раскраски графов на языке C++. В нем реализованы следующие алгоритмы:

### **Жадный алгоритм расскраски (seq_greedy)** 
Последовательный жадный алгоритм обращается ко всем вершинам в случайном порядке и присваивает каждой из них цвет в зависимости от цветов соседей. В частности, каждая вершина окрашивается наименьшим недостающим цветом из ее окрестности, пока не будут раскрашены все вершины.

В нашей реализации порядок доступа к вершинам рандомизируется при каждом запуске, что приводит к недетерминированной раскраске.

### **Алгоритм расскраски с выбором вершины максимальной степени (seq_ldf)** 
Последовательный алгоритм LDF обращается ко всем вершинам, начиная с вершин с наибольшей степенью, и присваивает каждой из них цвет так же, как это делает жадный алгоритм.

В нашей реализации вершины сортируются с помощью сортировки кучей по их степеням.

### **Рекурсивный алгоритм с выбором максимально независимого множества (rec_rlf)**
Алгоритм RLF раскрашивает вершины графа, создавая один цветовой класс за раз. Он определяет максимальный независимый набор вершин в графе, присваивает им одинаковый цвет, а затем удаляет эти вершины из графа. Этот процесс повторяется для оставшегося подграфа до тех пор, пока все вершины не будут раскрашены.

### **Параллельный алгоритм расскраски Джонса-Плассмана (par_jp)** 
Coming soon...

### **Параллельный алгоритм оптимизации раскраски по вершинам (color_op_ver)** 
Coming soon...

### **Параллельный алгоритм оптимизации раскраски по рёбрам (color_op_ed)** 
Coming soon...

Проект поддерживает параллельные вычисления и позволяет сохранять результаты в файл CSV. В нем также имеется скрипт на Python для визуализации результатов работы алгоритмов. Два из этих алгоритмов являются собственными и представляют собой уникальные методы раскраски графов. 