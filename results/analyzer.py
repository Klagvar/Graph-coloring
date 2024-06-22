import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # Считываем все результаты в один DataFrame
    df_list = []
    for entry in os.scandir('.'):
        if entry.path.endswith('.csv'):
            data = pd.read_csv(entry.path)
            df_list.append(data)
    df = pd.concat(df_list, ignore_index=True)

    # ПОСТРОЕНИЕ ГРАФИКОВ ДЛЯ ВСЕХ ГРАФОВ

    # Создаем фигуры для графиков
    fig1 = plt.figure(figsize=(15, 8))
    ax1 = fig1.add_subplot(111)
    fig2 = plt.figure(figsize=(15, 8))
    ax2 = fig2.add_subplot(111)

    # Для каждого метода раскраски
    for coloring_method in df.coloring_method.unique():
        # Группируем данные по имени графа
        df_vc = df[df.coloring_method == coloring_method].groupby(
            ['graph_name'])[['coloring_time', 'colors_used']].mean().reset_index()

        # Конвертируем имена графов в строки
        x = [str(i) for i in df_vc.graph_name.tolist()]
        y = df_vc.coloring_time.tolist()
        ax1.plot(x, y, '-o', label=coloring_method, linewidth=0.75)
        y = df_vc.colors_used.tolist()
        ax2.plot(x, y, '-o', label=coloring_method, linewidth=0.75)

    ax1.legend(prop={'size': 11})
    ax2.legend(prop={'size': 11})
    ax1.grid()
    ax2.grid()
    ax1.set_xlabel('Название графа', fontsize=13)
    ax2.set_xlabel('Название графа', fontsize=13)
    ax1.set_ylabel('Время выполнения в секундах', fontsize=13)
    ax2.set_ylabel('Количество использованных цветов', fontsize=13)
    ax1.set_title('Сравнение алгоритмов раскраски на 34 графах', fontsize=15)
    ax2.set_title('Сравнение алгоритмов раскраски на 34 графах', fontsize=15)
    ax1.set_yticks(range(0, int(max(ax1.get_yticks())-1)))
    ax1.set_xticks(range(len(x)))
    ax1.set_xticklabels(x, rotation=50, ha='right', fontsize=11)
    ax2.set_xticks(range(len(x)))
    ax2.set_xticklabels(x, rotation=50, ha='right', fontsize=11)
    ax1.tick_params(axis='x', labelsize=12)
    ax2.tick_params(axis='x', labelsize=12)
    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig('plots/all_time_comparison.png')
    fig2.savefig('plots/all_colors_comparison.png')

    # ПОСТРОЕНИЕ ГРАФИКОВ ДЛЯ ГРАФОВ RGG

    # Фильтруем DataFrame по графам RGG
    df_rgg = df[df.graph_name.str.contains('^rgg')]

    # Создаем фигуры для графиков по количеству вершин
    fig1 = plt.figure(figsize=(12, 8))
    ax1 = fig1.add_subplot(111)
    fig2 = plt.figure(figsize=(12, 8))
    ax2 = fig2.add_subplot(111)

    # Для каждого метода раскраски
    for coloring_method in df.coloring_method.unique():
        # Группируем данные по количеству вершин
        df_vc = df_rgg[df_rgg.coloring_method == coloring_method].groupby(
            ['vertex_count'])[['coloring_time', 'colors_used']].mean().reset_index()

        # Конвертируем количество вершин в строки
        x = [str(i) for i in df_vc.vertex_count.tolist()]
        y = df_vc.coloring_time.tolist()
        ax1.plot(x, y, '-o', label=coloring_method, linewidth=0.75)
        y = df_vc.colors_used.tolist()
        ax2.plot(x, y, '-o', label=coloring_method, linewidth=0.75)

    ax1.legend(prop={'size': 12})
    ax2.legend(prop={'size': 12})
    ax1.grid()
    ax2.grid()
    ax1.set_xlabel('Количество вершин', fontsize=13)
    ax2.set_xlabel('Количество вершин', fontsize=13)
    ax1.set_ylabel('Время выполнения в секундах', fontsize=13)
    ax2.set_ylabel('Количество использованных цветов', fontsize=13)
    ax1.set_title('Сравнение алгоритмов раскраски на случайных геометрических графах', fontsize=15)
    ax2.set_title('Сравнение алгоритмов раскраски на случайных геометрических графах', fontsize=15)
    ax1.set_yticks(np.arange(0, 420, 20))
    fig1.tight_layout()
    fig2.tight_layout()
    ax1.tick_params(axis='x', labelsize=12)
    ax2.tick_params(axis='x', labelsize=12)
    fig1.savefig('plots/rgg_time_comparison.png')
    fig2.savefig('plots/rgg_colors_comparison.png')

    # ПОСТРОЕНИЕ ГРАФИКОВ ПО КОЛИЧЕСТВУ РЁБЕР ДЛЯ ГРАФОВ RGG

    # Извлекаем количество рёбер из названия графа и добавляем его в DataFrame
    df_rgg['edge_count'] = df_rgg['graph_name'].str.extract(r'rgg_(\d+)\.')[0].astype(int)

    fig3 = plt.figure(figsize=(12, 8))
    ax3 = fig3.add_subplot(111)

    # Для каждого метода раскраски
    for coloring_method in df_rgg.coloring_method.unique():
        # Группируем данные по количеству рёбер
        df_ec = df_rgg[df_rgg.coloring_method == coloring_method].groupby(
            ['edge_count'])[['coloring_time']].mean().reset_index()

        # Конвертируем количество рёбер в строки
        x = [str(i) for i in df_ec.edge_count.tolist()]
        y = df_ec.coloring_time.tolist()
        ax3.plot(x, y, '-o', label=coloring_method, linewidth=0.75)

    ax3.legend(prop={'size': 12})
    ax3.grid()
    ax3.set_xlabel('Количество рёбер', fontsize=13)
    ax3.set_ylabel('Время выполнения в секундах', fontsize=13)
    ax3.set_title('Сравнение алгоритмов раскраски по количеству рёбер', fontsize=15)
    ax3.set_yticks(np.arange(0, 620, 20))
    ax3.set_xticks(range(len(x)))
    ax3.tick_params(axis='x', labelsize=12)
    fig3.tight_layout()
    fig3.savefig('plots/rgg_edge_comparison.png')

    # ПОСТРОЕНИЕ ГРАФИКА С СОРТИРОВКОЙ ПО КОЛИЧЕСТВУ ЦВЕТОВ ДЛЯ АЛГОРИТМА par_opt_ed

    # Создаем фигуру для графика
    fig4 = plt.figure(figsize=(15, 8))
    ax4 = fig4.add_subplot(111)

    # Группируем данные по имени графа для алгоритма par_opt_ed
    df_sorted = df[df.coloring_method == 'par_opt_ed'].groupby(
        ['graph_name'])[['colors_used']].mean().reset_index()
    
    # Сортируем данные по количеству использованных цветов
    df_sorted = df_sorted.sort_values(by='colors_used')
    sorted_graph_names = df_sorted.graph_name.tolist()

    # Для каждого метода раскраски
    for coloring_method in df.coloring_method.unique():
        df_vc = df[df.coloring_method == coloring_method]
        df_vc = df_vc[df_vc.graph_name.isin(sorted_graph_names)].groupby(
            ['graph_name'])[['colors_used']].mean().reindex(sorted_graph_names).reset_index()

        # Конвертируем имена графов в строки
        x = [str(i) for i in df_vc.graph_name.tolist()]
        y = df_vc.colors_used.tolist()
        ax4.plot(x, y, '-o', label=coloring_method, linewidth=0.75)

    ax4.legend(prop={'size': 11})
    ax4.grid()
    ax4.set_xlabel('Название графа', fontsize=13)
    ax4.set_ylabel('Количество использованных цветов', fontsize=13)
    ax4.set_title('Сравнение алгоритмов раскраски, отсортированное по количеству цветов в par_opt_ed', fontsize=15)
    ax4.set_xticks(range(len(x)))
    ax4.set_xticklabels(x, rotation=50, ha='right', fontsize=11)
    ax4.tick_params(axis='x', labelsize=12)
    fig4.tight_layout()
    fig4.savefig('plots/sorted_by_par_opt_ed_colors.png')
