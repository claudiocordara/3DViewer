#include "mainwindow.h"
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QtWidgets>
#include <Windows.h>
#include "scene3d.h"
#include "functions.h"

// The MainWindow class to visualization of 3D-objects and their processing control
MainWindow::MainWindow()
{
	// create the window controls
	// define the action object to menu items creation
	QAction * action;
	
	// create the 'File' menu which will provide the file loading/saving
	QMenu * menu = menuBar()->addMenu(tr("&File"));
	// create 'New Model' item
	menu->addAction(tr("New Model"), this, &MainWindow::openModel);
	// create 'Append Part' item
	menu->addAction(tr("Append Part"), this, &MainWindow::addPart);
	menu->addSeparator();
	// create 'Save Parts' item
	menu->addAction(tr("Save Parts"), this, &MainWindow::savePOFF);
	// create 'Save Single' item
	menu->addAction(tr("Save Single"), this, &MainWindow::saveSOFF);
	menu->addSeparator();
	// create 'Quit' item
	menu->addAction(tr("&Quit"), this, &QWidget::close);

	// create the 'Process' menu which will provide the start of different calculations
	menuActions = menuBar()->addMenu(tr("&Process"));
	// create 'Create Skeleton' item and set it as disabled
	action = menuActions->addAction(tr("Create Skeleton"), this, &MainWindow::Skeleton);
	action->setEnabled(false);
	// create 'Save separate parts' item and set it as disabled
	action = menuActions->addAction(tr("Save separate parts"), this, &MainWindow::sGroups);
	action->setEnabled(false);
	menuActions->addSeparator();
	// create 'Get new part' item and set it as disabled
	action = menuActions->addAction(tr("Get new part"), this, &MainWindow::nPart);
	action->setEnabled(false);
	menuActions->addSeparator();
	action = menuActions->addAction(tr("Split and Save by Skeleton"), this, &MainWindow::splitAndSaveBySkeleton);
	action->setEnabled(true);
	action = menuActions->addAction(tr("Split and Save by SDF"), this, &MainWindow::splitAndSaveBySDF);
	action->setEnabled(true);
	action = menuActions->addAction(tr("Split and Save by Skeleton and SDF"), this, &MainWindow::splitAndSaveBySkeletonAndSDF);
	action->setEnabled(true);


	mLayout = new QHBoxLayout();
	//QPushButton *button1 = new QPushButton("One");
	// create widget (Scene3D object) to show the 3D-objects
	widget = new Scene3D(this);
	mTree = new QTreeWidget();
	mTree->setColumnCount(1);

	//QList<QTreeWidgetItem *> items;
	//for (int i = 0; i < 10; ++i)
	//	items.append(new QTreeWidgetItem((QTreeWidget*)0, QStringList(QString("item: %1").arg(i))));
	//tree->insertTopLevelItems(0, items);
	//QTreeWidgetItem *itm1 = new QTreeWidgetItem(QStringList(QString("1")));
	//QTreeWidgetItem *itm2 = new QTreeWidgetItem(QStringList(QString("2")));
	//itm1->addChild(itm2);
	//tree->insertTopLevelItem(0, itm1);
	if (QTreeWidgetItem* header = mTree->headerItem()) {
		header->setText(0, "Segmentation");
	}
	else {
		mTree->setHeaderLabel("Segmentation");
	}

	QSizePolicy spLeft(QSizePolicy::Preferred, QSizePolicy::Preferred);
	spLeft.setHorizontalStretch(2);
	widget->setSizePolicy(spLeft);
	QSizePolicy spRight(QSizePolicy::Preferred, QSizePolicy::Preferred);
	spRight.setHorizontalStretch(1);
	mTree->setSizePolicy(spRight);


	mLayout->addWidget(widget);
	mLayout->addWidget(mTree);

	QWidget *central = new QWidget();
	// Set it as central widget of window
	setCentralWidget(central);
	centralWidget()->setLayout(mLayout);

	// create the 'Elements' menu which will provide the control of visibility of elements
	menuOptions = menuBar()->addMenu(tr("Elements"));
	// create the checker for axis
	action = menuOptions->addAction(tr("Axis"));
	action->setCheckable(true);
	action->setChecked(true);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions);
	// create the checker for solid's skeleton
	action = menuOptions->addAction(tr("Skeleton"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions);
	// create the checker for solid's wireframe
	action = menuOptions->addAction(tr("Wireframe"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions);
	// create the checker for solid's facets
	action = menuOptions->addAction(tr("Facets"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions);
	// create the checker for extracted segments
	action = menuOptions->addAction(tr("Parts"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions);

	// create the checker for colored segments
	action = menuOptions->addAction(tr("Colored segments"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions1);
	action = menuOptions->addAction(tr("Colored segments test"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions2);
	action = menuOptions->addAction(tr("Colored sdf"));
	action->setCheckable(true);
	action->setChecked(false);
	action->setEnabled(false);
	connect(action, &QAction::toggled, this, &MainWindow::setDockOptions3);


	menuTest = menuBar()->addMenu(tr("Test"));
	action = menuTest->addAction(tr("Segmentation By SDF"), this, &MainWindow::TestSegmentationBySDF);
	action = menuTest->addAction(tr("Segmentation By Skeleton and SDF"), this, &MainWindow::TestSegmentationBySkeletonAndSDF);
	action = menuTest->addAction(tr("Polyedra Decomposition"), this, &MainWindow::TestPolyedraDecomposition);

	mTree->clear();

}

// Open the STL and OFF files
void MainWindow::openModel()
{
	// define the CGAL mesh object
	Polyhedron mesh;
	// get the filename using the standard Qt file dialog
	QString fileName = QFileDialog::getOpenFileName(this, "Load", QDir::currentPath(), "Models (*.stl);;Models (*.off)");
	// extract the file extension
	QString ext = fileName.right(3);
	
	// if the file is STL
	if (QString::compare(ext, "stl", Qt::CaseInsensitive) == 0) {
		// check the type of STL (ascii or binary)
		int type = getStlFileFormat(fileName);
		char * offName;
		// if type is binary
		if (type == STL_BINARY)
			// convert binary STL to OFF
			offName = binStl2Off(fileName.toLatin1().data());
		// if type is ascii
		else if (type == STL_ASCII)
			// convert ascii STL to OFF
			offName = ascStl2Off(fileName.toLatin1().data());
		// unrecognized type of STL
		else {
			// show the error message
			QMessageBox msgBox;
			msgBox.setText("This file corrupt and cannot be loaded");
			msgBox.setInformativeText("STL Load");
			msgBox.setStandardButtons(QMessageBox::Ok);
			msgBox.exec();
			return;
		}
		// load OFF into mesh using CGAL
		std::ifstream input(offName);
		input >> mesh;
	}
	// if the file is OFF
	else {
		// load OFF into mesh using CGAL
		QByteArray ba = fileName.toLatin1();
		std::ifstream input(ba.data());
		input >> mesh;
	}
	// delete the existing widget
	delete widget;
	// create the widget with new mesh
	mLayout->removeWidget(widget);
	widget = new Scene3D(this);
	QSizePolicy spLeft(QSizePolicy::Preferred, QSizePolicy::Preferred);
	spLeft.setHorizontalStretch(2);
	widget->setSizePolicy(spLeft);
	widget->load(mesh);
	// set it as central widget of window
	mLayout->insertWidget(0, widget);
	//setCentralWidget(widget);
	setFigureOn();
}

// Skeleton creation start from GUI
void MainWindow::Skeleton()
{
	// call the method of Scene3D (see header)
	int mGroup = createSkeleton();
	// change the statement of 'Skeleton' and 'Colored segments' checker
	QList<QAction*> actions = menuOptions->actions();
	actions.at(1)->setChecked(true);
	actions.at(1)->setEnabled(true);
	actions.at(5)->setChecked(false);
	actions.at(5)->setEnabled(true);

	// disable the 'Create Skeleton' item
	actions = menuActions->actions();
	actions.at(0)->setEnabled(false);
	// check do we have extracted parts
	if (mGroup > 0)
		// enable 'Save separate parts' item
		actions.at(1)->setEnabled(true);
	// enable 'Get new part' item
	actions.at(3)->setEnabled(true);

	// refresh the 'elements visibility' variable
	setDockOptions();
}

// Save all the extracted segments
void MainWindow::sGroups()
{
	// call the method of Scene3D (see header)
	saveGroups();
}

// Try to extract the new segment
void MainWindow::nPart()
{
	// call the method of Scene3D
	widget->newPart();
	
	// use the 'Elements' menu
	QList<QAction*> actions = menuOptions->actions();
	// disable 'Skeleton' checker
	actions.at(1)->setChecked(false);
	actions.at(1)->setEnabled(false);
	
	// use the 'Process' menu
	actions = menuActions->actions();
	// enable 'New Model' item 
	actions.at(0)->setEnabled(true);
	// disable 'Append Part' item
	actions.at(1)->setEnabled(false);
	// disable 'Save Parts' item
	actions.at(3)->setEnabled(false);

	// refresh the 'elements visibility' variable
	setDockOptions();
}


void MainWindow::splitAndSaveBySkeleton() {
	widget->splitPartsBySkeleton();

	QList<QAction*> actions = menuOptions->actions();
	actions.at(5)->setEnabled(true);
}

void MainWindow::splitAndSaveBySDF() {
	widget->splitPartBySDF();

	QList<QAction*> actions = menuOptions->actions();
	actions.at(6)->setEnabled(true);
	actions.at(7)->setEnabled(true);
}


void MainWindow::splitAndSaveBySkeletonAndSDF() {
	widget->splitPartBySkeletonAndSDF();

	QList<QAction*> actions = menuOptions->actions();
	actions.at(6)->setEnabled(true);
	actions.at(7)->setEnabled(true);
}


// Save all #D-objects (main and extracted parts) to OFF files
void MainWindow::savePOFF()
{
	// get the filename using the standard Qt file dialog
	QString partName = QFileDialog::getSaveFileName(this, "Save", QDir::currentPath(), "Part (*.off)");
	// loop over the existing objects
	for (int part = 0; part < widget->tmesh.size(); part++) {
		// create the serial part of filename
		std::stringstream sc;
		sc.str(std::string());
		sc << std::setw(3) << std::setfill('0') << part << ".off";
		QString partN = QString::fromStdString(sc.str());
		// replace the extension with a serial part 
		QByteArray ba = partName.toLatin1();
		QString curName = QString(ba).replace(".off", partN);
		// call the method from functions module
		meshToOff(curName.toLatin1().data(), widget->tmesh[part]);
	}
}

// Save the single (main) 3D-object to OFF file
void MainWindow::saveSOFF()
{
	// get the filename using the standard Qt file dialog
	QString fileName = QFileDialog::getSaveFileName(this, "Save", QDir::currentPath(), "Part (*.off)");
	// call the method from functions module
	meshToOff(fileName.toLatin1().data(), mergePoly(widget->tmesh, 0, (int)widget->tmesh.size()));
}

// Add new part from file
void MainWindow::addPart()
{
	// define the CGAL mesh object
	Polyhedron mesh;
	// get the filename using the standard Qt file dialog
	QString fileName = QFileDialog::getOpenFileName(this, "Load", QDir::currentPath(), "Models (*.stl);;Models (*.off)");
	// extract the extension
	QString ext = fileName.right(3);

	// if the file is STL
	if (QString::compare(ext, "stl", Qt::CaseInsensitive) == 0) {
		// check the type of STL (ascii or binary)
		int type = getStlFileFormat(fileName);
		char * offName;
		// if type is binary
		if (type == STL_BINARY)
			// convert binary STL to OFF
			offName = binStl2Off(fileName.toLatin1().data());
		// if type is ascii
		else if (type == STL_ASCII)
			// convert ascii STL to OFF
			offName = ascStl2Off(fileName.toLatin1().data());
		// unrecognized type of STL
		else {
			// show the error message
			QMessageBox msgBox;
			msgBox.setText("This file corrupt and cannot be loaded");
			msgBox.setInformativeText("STL Load");
			msgBox.setStandardButtons(QMessageBox::Ok);
			msgBox.exec();
			return;
		}
		// load OFF into mesh using CGAL
		std::ifstream input(offName);
		input >> mesh;
	}
	// if the file is OFF
	else {
		// load OFF into mesh using CGAL
		QByteArray ba = fileName.toLatin1();
		std::ifstream input(ba.data());
		input >> mesh;
	}
	// add the loaded mesh as a new part to Scene3D object
	widget->add(mesh);
	// refresh the menu
	setFigureOn();
}

// Set the initial statement of menu items
void MainWindow::setFigureOn()
{
	// use the 'Elements' menu
	QList<QAction*> actions = menuOptions->actions();
	// enable and set on 'Wireframe' checker
	actions.at(2)->setChecked(true);
	actions.at(2)->setEnabled(true);
	// enable and set on 'Facets' checker
	actions.at(3)->setChecked(true);
	actions.at(3)->setEnabled(true);
	// enable and set on 'Parts' checker
	actions.at(4)->setChecked(true);
	actions.at(4)->setEnabled(true);

	// use the 'Process' menu
	actions = menuActions->actions();
	// enable 'Create Skeleton'
	actions.at(0)->setEnabled(true);
	// disable 'Save Parts'
	actions.at(1)->setEnabled(false);

	// refresh the 'elements visibility' variable
	setDockOptions();
}

// Create the 'elements visibility' variable according to checkers of 'Elements' menu
void MainWindow::setDockOptions()
{
	// use the 'Elements' menu
	QList<QAction*> actions = menuOptions->actions();
	// initiate variable with zero
	widget->showElem = 0;
	
	// set the mask of 'Axis' item
	if (actions.at(0)->isChecked())
		widget->showElem |= shAxis;
	// set the mask of 'Skeleton' item
	if (actions.at(1)->isChecked())
		widget->showElem |= shSkeleton;
	// set the mask of 'Wireframe' item
	if (actions.at(2)->isChecked())
		widget->showElem |= shWireframe;
	// set the mask of 'Facets' item
	if (actions.at(3)->isChecked())
		widget->showElem |= shFacets;
	// set the mask of 'Parts' item
	if (actions.at(4)->isChecked())
		widget->showElem |= shParts;
	// set the mask of 'Colored segments' item
	if (actions.at(5)->isChecked()) {
		widget->showElem |= shColors;
		widget->showElem &= ~segColors;
		widget->showElem &= ~sdfColors;
	}
	if (actions.at(6)->isChecked()) {
		int ntshColors = ~shColors;
		widget->showElem &= ~shColors;
		widget->showElem |= segColors;
		widget->showElem &= ~sdfColors;
	}
	if (actions.at(7)->isChecked()) {
		widget->showElem &= ~shColors;
		widget->showElem &= ~segColors;
		widget->showElem |= sdfColors;
	}

	// set the colors (to apply 'Colored segments' item)
	sColors();

	// update the showed elements
	widget->update();
}

void MainWindow::setDockOptions1() {
	QList<QAction*> actions = menuOptions->actions();
	if (actions.at(5)->isChecked()) {
		actions.at(6)->setChecked(false);
		actions.at(7)->setChecked(false);
		setDockOptions();
	}
}
void MainWindow::setDockOptions2() {
	QList<QAction*> actions = menuOptions->actions();
	if (actions.at(6)->isChecked()) {
		actions.at(5)->setChecked(false);
		actions.at(7)->setChecked(false);
		setDockOptions();
	}
}
void MainWindow::setDockOptions3() {
	QList<QAction*> actions = menuOptions->actions();
	if (actions.at(7)->isChecked()) {
		actions.at(5)->setChecked(false);
		actions.at(6)->setChecked(false);
		setDockOptions();
	}
}



int MainWindow::TestSegmentationBySDF() {
	int ret = widget->testSegmentationBySDF();
	if (ret == EXIT_SUCCESS) {
		QList<QAction*> actions = menuOptions->actions();
		actions.at(6)->setEnabled(true);
		actions.at(7)->setEnabled(true);
	}
	return ret;
}


int MainWindow::TestSegmentationBySkeletonAndSDF() {
	int ret = widget->testSegmentationBySkeletonAndSDF();
	if (ret == EXIT_SUCCESS) {
		QList<QAction*> actions = menuOptions->actions();
		actions.at(6)->setEnabled(true);
		actions.at(7)->setEnabled(true);
	}
	return ret;
}


int MainWindow::TestPolyedraDecomposition() {
	int ret = widget->testPolyedraDecomposition();
	return ret;
}

void MainWindow::clearTree() {
	mTree->clear();
	mItems.clear();
}


void MainWindow::addTreeItem(int _parentIndex, QString _text) {
	QTreeWidgetItem *itm = new QTreeWidgetItem(QStringList(_text));
	if (_parentIndex >= 0 && _parentIndex < mItems.size()) {
		mItems[_parentIndex]->addChild(itm);
	}
	else {
		mItems.append(itm);
		mTree->addTopLevelItem(itm);
	}
}
