# Homework4-5

**HW4**

Был создан метод visualise (self, fullness, output_file), принимающий на вход опцию отображения (делать вывод полным = anything - в узлах и на ребрах подписаны последовательности или сокращенным = short в узлах подписаны покрытия, на ребрах – длины и покрытия), а также файл, куда будет сохранен граф в формате dot.

Для k = 3 были сохранены изображения полного и сокращенного вывода в формате png (HW4_full_3.png, HW4_short_3.png).
Для k = 15 были сохранены полный и сокращенный выводы в формате dot (HW4_full_15.dot и HW4_short_15.dot).


**HW5**

Был создан метод compress, сжимающий исходный граф. После сжатия необходимо пересчитать покрытие граней графа встроенной функцией calculate_init_coverage.

В качестве проверки работоспособности, сначала были построены и сохранены в формате png сжатые и не сжатые графы для последовательности ACTACGTACGTACG$:

k = 3

my_graph = Graph(k)

my_graph.add_read('ACTACGTACGTACG$')

my_graph.calc_init_edge_coverage()

my_graph.visualize("full", "/test_full.png")

my_graph.compress()

my_graph.visualize("full", "test_full_compressed.png")

Рисунки прикреплены в директории. 

Кроме того, была проверена работоспособность на первых ста ридах из предложенного файла hw_4_5_dataset.fasta для k = 15 (код в тексте скрипта). Полученные аутпуты были сохранены и приложены в репозиторий как для не сжатых графов, так и для сжатых в полном и укороченном формате - файлы названы hw_5_15full.dot, hw_5_15short.dot для не сжатых графов и hw_5_compressed_15full.dot, hw_5_comressed_15short.dot для сжатых графов. 


