# imports extern

# imports intern
import src.argument_parser as pars
import src.output_browser as gui
import src.output_terminal as pix


def launch_me():
    if pars.args.display:
        # launch dash version
        gui.run_gui()
    else:
        # launch terminal as default
        pix.use_terminal()


if __name__ == "__main__":
    launch_me()
