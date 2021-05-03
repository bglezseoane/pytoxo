# -*- coding: utf-8 -*-

###########################################################
# PyToxo
#
# A Python library for calculating penetrance tables of any
# bivariate epistasis model.
#
# Copyright 2021 Borja González Seoane
#
# Contact: borja.gseoane@udc.es
###########################################################

"""Graphical user interface entry point."""

import PySimpleGUI as sg


def main():
    # Style
    sg.theme("DarkBlue4")
    window_general_font = "Verdana, 20"

    # Model table components
    headings = [
        "Genotype",
        "Penetrance expression",
        "Calculated penetrance probability",
    ]
    header = [[sg.Text(h, size=(20, 2), pad=(0, 0)) for h in headings]]

    input_rows = [
        [sg.Input(size=(20, 1), pad=(0, 0)) for col in range(len(headings))]
        for row in range(10)
    ]

    # Main layout
    layout = [
        [
            sg.Combo(
                ("Heritability", "Prevalence"),
                key="prev_or_her_cb",
                size=(10, 1),
                readonly=True,
            ),
            sg.InputText(key="prev_or_her_in", size=(5, 1)),
        ],
        [
            sg.Text("MAFs (separated with commas)"),
            sg.InputText(key="mafs_in", size=(5, 1)),
        ],
        [
            sg.Button("Calculate table"),
        ],
        [header + input_rows],
    ]

    window = sg.Window(
        "PyToxo GUI",
        layout,
        default_element_size=(40, 1),
        grab_anywhere=False,
        font=window_general_font,
    )

    while True:
        event, values = window.read()

        if event == "SaveSettings":
            filename = sg.popup_get_file("Save Settings", save_as=True, no_window=True)
            window.SaveToDisk(filename)
            # save(values)
        elif event == "LoadSettings":
            filename = sg.popup_get_file("Load Settings", no_window=True)
            window.LoadFromDisk(filename)
            # load(form)
        elif event in ("Exit", None):
            break

    window.close()


if __name__ == "__main__":
    main()
