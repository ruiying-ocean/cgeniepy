    def skill_score(self, cost_function="m_score", table_styler=True, *args, **kwargs):
        "summarised model skill score compared to modern observations"

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())
        df = {
            "Biomass": [
                ForamType(i, self.model_path)
                .biomass_c()
                ._run_method(method=cost_function, *args, **kwargs)
                for i in foram_abbrev
            ],
            "Carbon Export": [
                ForamType(i, self.model_path)
                .export_c()
                ._run_method(method=cost_function, *args, **kwargs)
                for i in foram_abbrev
            ],
            "Relative Abundance": [
                ForamType(i, self.model_path)
                .export_c()
                .proportion()
                ._run_method(method=cost_function, *args, **kwargs)
                for i in foram_abbrev
            ],
        }

        df = DataFrame(df, index=foram_fullname)

        df["Column Total"] = df.sum(axis=1)
        df.loc["Row Total", :] = df.sum(axis=0)

        if table_styler:
            df = df.style.set_caption(
                f"{cost_function} across foraminifer groups and variables compared to modern observation"
            ).text_gradient(cmap="viridis", subset=(df.index[0:4], df.columns[0:3]))

        return df

    def summarise(
        self,
        stat="nanmean",
        diff=False,
        diff_method="percentage",
        table_styler=True,
        *args,
        **kwargs,
    ):
        """
        summarise basic statistics of foraminifer groups
        """

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())

        dic = {
            "Biomass": [
                ForamType(i, self.model_path).biomass_c()._run_method(method=stat)
                for i in foram_abbrev
            ],
            "Carbon Export": [
                ForamType(i, self.model_path).export_c()._run_method(method=stat)
                for i in foram_abbrev
            ],
            "Relative Abundance": [
                ForamType(i, self.model_path)
                .export_c()
                .proportion()
                ._run_method(method=stat)
                for i in foram_abbrev
            ],
        }
        df = DataFrame(dic, index=foram_fullname)

        df["Column Total"] = df.sum(axis=1)
        df.loc["Row Total", :] = df.sum(axis=0)

        if diff:
            obs = obs_stat_bytype(type=stat, *args, **kwargs)
            diff = df - obs

            if diff_method == "percentage":
                diff_perc = diff / obs
                return diff_perc
            elif diff_method == "absolute":
                return diff

        if table_styler:
            df = df.style.set_caption(
                f"{stat} across foraminifer groups"
            ).text_gradient(cmap="viridis")

        return df

    def foram_ldg(self, legend=True):
        "plot ldg: latitudinal diversity gradient"
        ## TODO: change colormap; change line width; add minor ticks
        ## add background grids
        foram_array_list = []
        foram_abbrev = list(foram_names().keys())
        foram_fullname = list(foram_names().values())

        for foram in foram_abbrev:
            array_1d = self.select_foram(foram).export_c().proportion().mean(axis=1)
            foram_array_list.append(array_1d)

        lat = GENIE_lat()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)

        for n in range(4):
            ax.plot(lat, foram_array_list[n], label=foram_fullname[n])

        if legend:
            ax.legend(
                loc="lower center",
                ncol=2,
                bbox_to_anchor=(0.5, -0.35),
                edgecolor="black",
            )

        ax.set_xlabel("Latitude")
        ax.set_ylabel("Relative abundance")
        ax.tick_params(axis="y", direction="in")
        ax.tick_params(axis="x", direction="in")

        return ax

    def plot_currents(self, z_level=0):
        """
        plot velocity field of horizontal currents, default is the toppest layer
        Not able to add upon other field map yet.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

        v = self.get_var("phys_v").normalise_longitude().isel(zt=z_level)
        u = self.get_var("phys_u").normalise_longitude().isel(zt=z_level)
        m = np.hypot(u, v)

        lat = v.lat
        lon = v.lon
        lon2d, lat2d = np.meshgrid(lon, lat)

        ax.set_global()
        ax.stock_img()
        ax.coastlines()
        ax.quiver(lon2d, lat2d, m, u, v, transform=ccrs.PlateCarree())

        return ax

    def barplot_comparison(self, *args, **kwargs) -> plt.axes:

        """
        Overview barplot of biomass and POC export compared to observed data
        :returns: matplotlib axes object
        """

        fname = list(foram_names().values())
        foram_abbrev = foram_names().keys()

        model_biomass_mean = [
            self.select_foram(i).biomass_c().nanmean() for i in foram_abbrev
        ]
        model_export_mean = [
            self.select_foram(i).export_c().nanmean() for i in foram_abbrev
        ]
        model_biomass_se = [self.select_foram(i).biomass_c().se() for i in foram_abbrev]
        model_export_se = [self.select_foram(i).export_c().se() for i in foram_abbrev]
        model_export_sum = [
            self.select_foram(i).export_c().sum().magnitude for i in foram_abbrev
        ]
        model_biomass_sum = [
            self.select_foram(i).biomass_c().sum().magnitude for i in foram_abbrev
        ]
        obs_biomass_mean = obs_stat_bysource("tow", *args, **kwargs).loc[:, "mean"]
        obs_biomass_se = obs_stat_bysource("tow", *args, **kwargs).loc[:, "se"]
        obs_export_mean = obs_stat_bysource("trap", *args, **kwargs).loc[:, "mean"]
        obs_export_se = obs_stat_bysource("trap", *args, **kwargs).loc[:, "se"]

        model_biomass_sum = DataFrame({"group": fname, "value": model_biomass_sum})
        model_export_sum = DataFrame({"group": fname, "value": model_export_sum})

        data_to_plot = [
            [
                [
                    model_biomass_mean,
                    model_biomass_se,
                    obs_biomass_mean,
                    obs_biomass_se,
                ],
                [model_export_mean, model_export_se, obs_export_mean, obs_export_se],
            ],
            [model_biomass_sum, model_export_sum],
        ]
        ## Plot starts
        fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)
        bar_width = 0.3

        # The x position of bars
        x = np.arange(4)
        xlabels = [w.replace(" ", "\n") for w in fname]

        for i in range(2):
            for j in range(2):
                axes[i, j].yaxis.set_minor_locator(AutoMinorLocator(4))
                axes[i, j].set_axisbelow(True)
                axes[i, j].yaxis.grid(color="gray", linestyle="dashed")
                if i == 0:
                    axes[i, j].bar(
                        x - bar_width / 2,
                        data_to_plot[i][j][0],
                        width=bar_width,
                        color=sns.color_palette("Set1")[0],
                        edgecolor="black",
                        yerr=data_to_plot[i][j][1],
                        capsize=7,
                        label="model",
                    )

                    axes[i, j].bar(
                        x + bar_width / 2,
                        data_to_plot[i][j][2],
                        width=bar_width,
                        color=sns.color_palette("Set1")[1],
                        edgecolor="black",
                        yerr=data_to_plot[i][j][3],
                        capsize=7,
                        label="obs",
                    )

                    axes[i, j].set_xticks(x)
                    axes[i, j].set_xticklabels(xlabels, rotation=45, ha="right")
                    axes[i, j].legend()
                else:
                    sns.barplot(
                        data=data_to_plot[i][j],
                        x="group",
                        y="value",
                        ax=axes[i, j],
                        edgecolor="black",
                        palette="deep",
                    )
                    axes[i, j].set_xticklabels(xlabels, rotation=45, ha="right")
                    axes[i, j].set_xlabel("")
                    axes[i, j].set_xticks(x)
                    set_sns_barwidth(axes[i, j], bar_width)

        axes[0, 0].set_ylabel(r"mmol C m$^{-3}$")
        axes[0, 1].set_ylabel(r"mmol C m$^{-2}$ d$^{-1}$")

        axes[0, 0].set_title("(a)    global biomass mean/se", loc="left")
        axes[0, 1].set_title("(b)    global POC export mean/se", loc="left")
        axes[1, 0].set_title("(c)    globally integrated biomass", loc="left")
        axes[1, 1].set_title("(d)    globally integrated EP", loc="left")

        axes[1, 0].set_ylabel("Gt C")
        axes[1, 1].set_ylabel(r"Gt C yr$^{-1}$")

        fig.tight_layout()

        return axes


def obs_stat_bysource(source, *args, **kwargs) -> pd.DataFrame:
    obs = []
    foram_abbrev = foram_names().keys()
    foram_fullnames = foram_names().values()[1]

    for i in foram_abbrev:
        tmp = obs_data(source=source, var=i, stat="Yes", *args, **kwargs)
        obs.append(tmp)

    table = pd.DataFrame(obs, index=foram_fullnames, columns=["mean", "sd", "se"])
    return table


def obs_stat_bytype(type, *args, **kwargs) -> pd.DataFrame:
    tow = obs_stat_bysource("net", *args, **kwargs).loc[:, type]
    trap = obs_stat_bysource("trap", *args, **kwargs).loc[:, type]
    core = obs_stat_bysource("coretop", *args, **kwargs).loc[:, type]
    # combination
    data = pd.concat([tow, trap, core], axis=1)
    data.columns = [
        "Biomass(mmol C/m3)",
        "Carbon Export (mmol C/m3/d)",
        "Relative Abundance",
    ]

    return data
